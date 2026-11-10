"""Physically meaningful diagnostics of a wind-driven gyre.

All integrals are discrete sums over interior nodes times ``dx dy``, the
same inner product under which the five-point Laplacian is symmetric. With
that choice the discrete energy and vorticity budgets close to round-off
(no quadrature error), so their residuals measure how far a state is from
steady rather than discretization error.

Quantities
----------
Energetics
    ``KE = -1/2 sum(psi zeta) dA = 1/2 int |grad psi|^2`` and enstrophy
    ``Z = 1/2 sum(zeta^2) dA``.
Energy budget
    ``dKE/dt = W + D_drag + D_visc + A_adv + A_beta`` with wind work
    ``W = -sum(psi F)``, drag ``D_drag = dS sum(psi zeta) = -2 dS KE``,
    viscous dissipation ``D_visc = -dM^3 sum(psi lap zeta)``, and advective
    and beta terms that vanish identically for the Arakawa Jacobian and the
    centered ``d/dx`` (kept to show that they do).
Vorticity budget
    Integrated over the basin, over the latitude band of the diagnosed gyre,
    and over the western strip and the remaining interior of that band: wind ``sum F``, beta ``-b sum psi_x``, drag ``-dS sum zeta``,
    viscosity ``dM^3 sum lap zeta`` and advection ``-dI^2 sum J``. The
    basin-integrated beta term is zero, so the wind's vorticity input can
    only be removed by friction, and the strip budgets show where that
    happens (the western boundary layer).
Boundary current
    At the latitude of the gyre center: peak speed and position, e-folding
    width and counter-current position (compared with boundary-layer
    theory), transport of the boundary current and its ratio to the
    Sverdrup prediction ``psi_I`` at the outer edge of the current. Mass
    conservation makes the current's transport equal to the southward
    interior transport east of that edge, so the ratio measures how closely
    the interior follows Sverdrup: below one when interior friction is not
    negligible, above one when a recirculation adds transport.
Resolution
    Grid points per boundary-layer decay scale, and the grid-scale fraction
    of the vorticity (norm of the part removed by a 3 x 3 running mean over
    the norm of the field). Values above about 0.3 mean that the smallest
    resolved scales are contaminated by two-grid-interval noise.
Asymmetry and concentration
    Intensification index (peak western speed over peak eastern speed),
    zonal kinetic-energy centroid (0.5 for a symmetric gyre), and the
    normalized Shannon entropy of the kinetic-energy distribution over the
    basin and of its zonal marginal. Western intensification concentrates
    energy into a narrow strip and lowers both entropies.
"""

from __future__ import annotations

import math

import numpy as np

from .analytic import measure_widths, sverdrup_interior, theory_widths
from .forcing import WindForcing
from .grid import Grid
from .operators import arakawa_jacobian, ddx, laplacian, velocity
from .params import Params


def _sum(f: np.ndarray, grid: Grid) -> float:
    return float(np.sum(f) * grid.cell_area)


def kinetic_energy(psi: np.ndarray, zeta: np.ndarray, grid: Grid) -> float:
    """``KE = -1/2 sum_interior(psi zeta) dx dy`` (equals ``1/2 int |u|^2`` discretely)."""
    return -0.5 * _sum(psi[1:-1, 1:-1] * zeta[1:-1, 1:-1], grid)


def enstrophy(zeta: np.ndarray, grid: Grid) -> float:
    return 0.5 * _sum(zeta[1:-1, 1:-1] ** 2, grid)


def normalized_entropy(weights: np.ndarray) -> float:
    """Shannon entropy of a non-negative weight array divided by ``log N``.

    Returns 1 for a uniform distribution and tends to 0 as all weight
    collects in a single cell.
    """
    w = np.asarray(weights, float).ravel()
    w = w[w > 0]
    if w.size <= 1:
        return 0.0
    q = w / w.sum()
    return float(-(q * np.log(q)).sum() / math.log(np.asarray(weights).size))


def gyre_latitude_index(psi: np.ndarray) -> int:
    """Index ``j`` of the row containing the largest ``|psi|`` (gyre center)."""
    return int(np.unravel_index(np.argmax(np.abs(psi)), psi.shape)[1])


def gyre_band(psi: np.ndarray, j: int) -> tuple[int, int]:
    """Latitude band ``[j0, j1]`` (node indices) of the gyre that contains row ``j``.

    The band is the contiguous set of rows on which the streamfunction at the
    longitude of the gyre center keeps the sign it has at row ``j``; for a
    double gyre it stops at the inter-gyre boundary.
    """
    i = int(np.argmax(np.abs(psi[:, j])))
    col = psi[i, :]
    s = np.sign(col[j])
    j0 = j
    while j0 > 0 and np.sign(col[j0 - 1]) == s:
        j0 -= 1
    j1 = j
    while j1 < len(col) - 1 and np.sign(col[j1 + 1]) == s:
        j1 += 1
    return j0, j1


def tendency_terms(psi: np.ndarray, zeta: np.ndarray, grid: Grid, p: Params,
                   F: np.ndarray) -> dict:
    """Interior arrays of every term on the right-hand side of the vorticity equation."""
    terms = {
        "wind": F[1:-1, 1:-1],
        "beta": -p.beta_factor * ddx(psi, grid.dx),
        "drag": -p.delta_s * zeta[1:-1, 1:-1],
        "viscosity": p.dm3 * laplacian(zeta, grid.dx, grid.dy),
        "advection": -p.di2 * arakawa_jacobian(psi, zeta, grid.dx, grid.dy),
    }
    return terms


def compute(psi: np.ndarray, zeta: np.ndarray, grid: Grid, p: Params, forcing: WindForcing,
            t: float = 0.0, strip: float | None = None, full: bool = True) -> dict:
    """Compute the full diagnostic set for one state.

    Parameters
    ----------
    psi, zeta:
        Full-grid streamfunction and vorticity (walls included).
    strip:
        Width of the western strip used in the regional vorticity budget.
        Default: three boundary-layer decay scales, at least four grid
        cells and at most half the basin.
    full:
        If False, skip the more expensive regional budgets and width fits.
    """
    F = forcing.curl_at(grid.X, grid.Y, t)
    d: dict = {"t": float(t)}
    ke = kinetic_energy(psi, zeta, grid)
    d["KE"] = ke
    d["enstrophy"] = enstrophy(zeta, grid)
    d["psi_max"] = float(psi.max())
    d["psi_min"] = float(psi.min())

    terms = tendency_terms(psi, zeta, grid, p, F)
    pi = psi[1:-1, 1:-1]
    d["energy_wind"] = -_sum(pi * terms["wind"], grid)
    d["energy_drag"] = -_sum(pi * terms["drag"], grid)
    d["energy_viscous"] = -_sum(pi * terms["viscosity"], grid)
    d["energy_advection"] = -_sum(pi * terms["advection"], grid)
    d["energy_beta"] = -_sum(pi * terms["beta"], grid)
    d["dKE_dt"] = (d["energy_wind"] + d["energy_drag"] + d["energy_viscous"]
                   + d["energy_advection"] + d["energy_beta"])
    scale = max(abs(d["energy_wind"]), 1e-300)
    d["energy_imbalance"] = d["dKE_dt"] / scale

    u, v = velocity(psi, grid.dx, grid.dy)
    e = 0.5 * (u**2 + v**2)
    d["max_speed"] = float(np.sqrt(2 * e.max()))
    d["ke_entropy"] = normalized_entropy(e)
    ex = e.sum(axis=1)
    d["ke_zonal_entropy"] = normalized_entropy(ex)
    d["ke_centroid_x"] = float((grid.x * ex).sum() / max(ex.sum(), 1e-300) / grid.Lx)

    j = gyre_latitude_index(psi)
    d["section_y"] = float(grid.y[j])
    row_v = v[:, j]
    w = measure_widths(grid.x, row_v)
    d["wbc_v_peak"] = w["v_peak"]
    d["wbc_x_peak"] = w["x_peak"]
    d["wbc_e_fold"] = w["e_fold"]
    d["wbc_zero_crossing"] = w["zero_crossing"]
    east = grid.x >= 0.5 * grid.Lx
    v_east = float(np.max(np.abs(row_v[east]))) if np.any(east) else math.nan
    d["intensification_index"] = abs(w["v_peak"]) / v_east if v_east > 0 else math.inf
    # boundary-current transport: psi where v first changes sign (edge of the current)
    sgn = np.sign(w["v_peak"]) if w["v_peak"] != 0 else 1.0
    half = grid.x <= 0.5 * grid.Lx
    row_psi = psi[:, j]
    if math.isfinite(w["zero_crossing"]):
        tr = float(np.interp(w["zero_crossing"], grid.x, row_psi))
    else:
        tr = float(row_psi[half][np.argmax(sgn * row_psi[half])])
    d["wbc_transport"] = tr
    if p.beta_factor > 0:
        x_edge = w["zero_crossing"] if math.isfinite(w["zero_crossing"]) else 0.0
        d["wbc_edge_x"] = float(x_edge)
        try:
            sv = float(sverdrup_interior(np.array(x_edge), np.array(grid.y[j]), p, forcing)) \
                * forcing.modulation(t)
        except ValueError:
            sv = math.nan
        d["sverdrup_transport"] = sv
        # the ratio is undefined where psi_I nearly vanishes (edge near the eastern
        # wall or a latitude of zero wind curl, common in eddying states)
        ref = abs(float(sverdrup_interior(np.array(0.0), np.array(grid.y[j]), p, forcing)))
        ok = math.isfinite(sv) and abs(sv) > 0.05 * max(ref, 1e-300)
        d["transport_ratio"] = tr / sv if ok else math.nan
    else:
        d["sverdrup_transport"] = math.nan
        d["transport_ratio"] = math.nan

    if not full:
        return d

    # grid-scale content of the vorticity: an under-resolved viscous scale shows up as
    # two-grid-interval noise that a 3x3 running mean removes
    zi = zeta[2:-2, 2:-2]
    nz = float(np.linalg.norm(zi))
    if nz > 0:
        from scipy.ndimage import uniform_filter
        d["grid_scale_fraction"] = float(np.linalg.norm(zi - uniform_filter(zi, 3)) / nz)
    else:
        d["grid_scale_fraction"] = 0.0

    # theory widths and resolution
    if p.beta_factor > 0 and p.has_friction:
        tw = theory_widths(p)
        d["theory_e_fold"] = tw["e_fold"]
        d["theory_zero_crossing"] = tw["zero_crossing"]
        d["theory_decay_scale"] = tw["decay_scale"]
        d["points_per_layer"] = tw["decay_scale"] / grid.dx
    else:
        d["theory_e_fold"] = d["theory_zero_crossing"] = d["theory_decay_scale"] = math.nan
        d["points_per_layer"] = math.nan

    # regional vorticity budgets
    if strip is None:
        base = d["theory_decay_scale"] if math.isfinite(d["theory_decay_scale"]) else 0.05
        strip = min(max(3.0 * base, 4 * grid.dx), 0.5 * grid.Lx)
    xi = grid.x[1:-1]
    j0, j1 = gyre_band(psi, j)
    band = np.zeros(grid.ny - 1, bool)
    band[max(j0 - 1, 0):j1] = True            # interior row index = j - 1
    west_x = (xi <= strip)[:, None]
    in_band = band[None, :]
    regions = {
        "basin": np.ones((grid.nx - 1, grid.ny - 1), bool),
        "gyre": np.broadcast_to(in_band, (grid.nx - 1, grid.ny - 1)),
        "west": west_x & in_band,
        "interior": (~west_x) & in_band,
    }
    budget = {}
    for rname, mask in regions.items():
        budget[rname] = {k: _sum(np.where(mask, arr, 0.0), grid) for k, arr in terms.items()}
        budget[rname]["tendency"] = sum(budget[rname].values())
    d["strip_width"] = float(strip)
    d["gyre_band"] = (float(grid.y[j0]), float(grid.y[j1]))
    d["vorticity_budget"] = budget
    d["regime"] = p.regime()
    return d


def describe(d: dict, scales=None) -> str:
    """Human-readable summary of :func:`compute` output.

    If ``scales`` (:class:`~dewi_lanjar.params.DimensionalScales`) is given,
    transports are also shown in sverdrups, widths in km and speeds in m/s.
    """
    def dim(val, kind):
        if scales is None or not math.isfinite(val):
            return ""
        if kind == "Sv":
            return f"  ({scales.sverdrups(val):.3g} Sv)"
        if kind == "km":
            return f"  ({scales.km(val):.3g} km)"
        if kind == "m/s":
            return f"  ({scales.m_per_s(val):.3g} m/s)"
        if kind == "days":
            return f"  ({scales.days(val):.3g} days)"
        return ""

    L = []
    if "regime" in d:
        L.append(f"Regime: {d['regime']['label']}")
    L.append(f"time                      t = {d['t']:.4g}{dim(d['t'], 'days')}")
    L.append(f"kinetic energy           KE = {d['KE']:.6g}")
    L.append(f"enstrophy                 Z = {d['enstrophy']:.6g}")
    L.append(f"streamfunction range        [{d['psi_min']:.4g}, {d['psi_max']:.4g}]"
             f"{dim(max(abs(d['psi_min']), abs(d['psi_max'])), 'Sv')}")
    L.append("")
    L.append(f"Western boundary current at y = {d['section_y']:.3f}")
    L.append(f"  peak velocity            v = {d['wbc_v_peak']:.4g}{dim(d['wbc_v_peak'], 'm/s')}")
    L.append(f"  peak position            x = {d['wbc_x_peak']:.4g}{dim(d['wbc_x_peak'], 'km')}")
    L.append(f"  e-folding width            {d['wbc_e_fold']:.4g}{dim(d['wbc_e_fold'], 'km')}"
             + (f"   theory {d['theory_e_fold']:.4g}" if math.isfinite(d.get('theory_e_fold', math.nan)) else ""))
    L.append(f"  counter-current at         {d['wbc_zero_crossing']:.4g}"
             + (f"   theory {d['theory_zero_crossing']:.4g}"
                if math.isfinite(d.get('theory_zero_crossing', math.nan)) else ""))
    L.append(f"  transport                  {d['wbc_transport']:.4g}{dim(d['wbc_transport'], 'Sv')}")
    L.append(f"  Sverdrup psi_I at edge     {d['sverdrup_transport']:.4g}")
    L.append(f"  transport ratio            {d['transport_ratio']:.4g}   (1 = Sverdrup closure)")
    L.append("")
    L.append("Asymmetry and concentration")
    L.append(f"  intensification index      {d['intensification_index']:.4g}   (1 = symmetric)")
    L.append(f"  KE zonal centroid x/Lx     {d['ke_centroid_x']:.4g}   (0.5 = symmetric)")
    L.append(f"  KE entropy (basin)         {d['ke_entropy']:.4g}   (1 = uniform)")
    L.append(f"  KE entropy (zonal)         {d['ke_zonal_entropy']:.4g}")
    L.append("")
    L.append("Energy budget  dKE/dt = W + D_drag + D_visc + A_adv + A_beta")
    for k in ("energy_wind", "energy_drag", "energy_viscous", "energy_advection", "energy_beta"):
        L.append(f"  {k:<18s} {d[k]: .6e}")
    L.append(f"  {'dKE/dt':<18s} {d['dKE_dt']: .6e}   (relative {d['energy_imbalance']:.2e})")
    if "vorticity_budget" in d:
        L.append("")
        y0, y1 = d["gyre_band"]
        L.append(f"Vorticity budget: gyre band y in [{y0:.3g}, {y1:.3g}], "
                 f"western strip x < {d['strip_width']:.3g}")
        cols = ("basin", "gyre", "west", "interior")
        L.append(f"  {'term':<11s}" + "".join(f"{r:>13s}" for r in cols))
        for term in ("wind", "beta", "drag", "viscosity", "advection", "tendency"):
            L.append(f"  {term:<11s}" + "".join(
                f"{d['vorticity_budget'][r][term]:>13.4e}" for r in cols))
    if "grid_scale_fraction" in d:
        L.append("")
        g_ = d["grid_scale_fraction"]
        flag = "" if g_ < 0.3 else "   WARNING: grid-scale noise, refine the grid or raise delta_m"
        L.append(f"Grid-scale vorticity fraction: {g_:.3f}{flag}")
    if math.isfinite(d.get("points_per_layer", math.nan)):
        L.append("")
        flag = "" if d["points_per_layer"] >= 3 else "   WARNING: boundary layer under-resolved"
        L.append(f"Resolution: {d['points_per_layer']:.3g} grid points per boundary-layer decay scale{flag}")
    return "\n".join(L)
