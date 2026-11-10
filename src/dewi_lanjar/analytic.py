"""Exact and asymptotic solutions of the steady linear problem.

The steady linear vorticity equation

.. math::

    \\delta_M^3 \\nabla^4\\psi - \\delta_S \\nabla^2\\psi - b\\,\\psi_x = -F(y)

with ``psi = 0`` on the walls and either free-slip (``zeta = 0``) or
no-slip (``psi_n = 0``) conditions separates exactly when the curl is a sum
of sines ``F = sum a_n sin(k_n y)``, ``k_n = n pi / Ly``. Writing
``psi = sum X_n(x) sin(k_n y)`` gives for each mode the constant-coefficient
ODE

.. math::

    \\delta_M^3 (D^2 - k^2)^2 X - \\delta_S (D^2 - k^2) X - b D X = -a ,

whose general solution is a constant particular part plus up to four
exponentials ``exp(m x)``. The roots ``m`` solve the characteristic quartic
(quadratic when ``delta_M = 0``). The constants follow from the wall
conditions. :func:`exact_linear` evaluates this solution to machine
precision on any set of points; it is the reference against which the
finite-difference solvers are verified.

The module also provides the Sverdrup interior and the classical
boundary-layer profiles of Stommel (1948) and Munk (1950).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .forcing import WindForcing
from .params import Params, boundary_layer_roots


# ---------------------------------------------------------------------------
# exact separable solution
# ---------------------------------------------------------------------------
@dataclass
class ModeSolution:
    """Exact solution ``X(x)`` of one sine mode.

    ``X(x) = Xp + sum_r c_r exp(m_r (x - x_r))`` with ``x_r = Lx`` for roots
    with positive real part and ``x_r = 0`` otherwise, so that no exponential
    exceeds one inside the basin (no overflow, good conditioning).
    """

    k: float
    a: float
    Xp: float
    roots: np.ndarray
    anchors: np.ndarray
    coeffs: np.ndarray

    def derivative(self, x: np.ndarray, order: int = 0) -> np.ndarray:
        x = np.asarray(x, float)
        e = np.exp(np.multiply.outer(x, self.roots) - self.roots * self.anchors)
        val = (e * (self.coeffs * self.roots**order)).sum(axis=-1).real
        if order == 0:
            val = val + self.Xp
        return val


def characteristic_roots(k: float, p: Params) -> np.ndarray:
    """Roots ``m`` of the characteristic polynomial for wavenumber ``k``.

    Quartic ``dM^3 m^4 - (2 dM^3 k^2 + dS) m^2 - b m + (dM^3 k^4 + dS k^2) = 0``,
    or the quadratic ``dS m^2 + b m - dS k^2 = 0`` when ``dM = 0``.
    """
    dm3, ds, b = p.dm3, p.delta_s, p.beta_factor
    if dm3 > 0:
        coeffs = [dm3, 0.0, -(2 * dm3 * k**2 + ds), -b, dm3 * k**4 + ds * k**2]
    elif ds > 0:
        coeffs = [ds, b, -ds * k**2]
    else:
        raise ValueError("the steady problem needs friction: delta_s > 0 or delta_m > 0")
    return np.roots(coeffs).astype(complex)


def solve_mode(k: float, a: float, p: Params) -> ModeSolution:
    """Exact solution of one sine mode of the steady linear problem."""
    roots = characteristic_roots(k, p)
    scale = max(1.0, float(np.max(np.abs(roots))))
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            if abs(roots[i] - roots[j]) < 1e-9 * scale:
                raise ValueError(
                    "repeated characteristic root (degenerate parameters); "
                    "use the finite-difference solver instead"
                )
    Lx = p.Lx
    Xp = -a / (p.dm3 * k**4 + p.delta_s * k**2)
    anchors = np.where(roots.real > 0, Lx, 0.0)

    conds = [(0.0, 0, -Xp), (Lx, 0, -Xp)]
    if p.dm3 > 0:
        d = 2 if p.bc == "free-slip" else 1
        conds += [(0.0, d, 0.0), (Lx, d, 0.0)]
    A = np.empty((len(conds), len(roots)), complex)
    rhs = np.empty(len(conds), complex)
    for row, (xb, d, val) in enumerate(conds):
        A[row] = roots**d * np.exp(roots * (xb - anchors))
        rhs[row] = val
    # column scaling keeps derivative rows (m^2 ~ delta^-2) well conditioned
    colscale = np.max(np.abs(A), axis=0)
    coeffs = np.linalg.solve(A / colscale, rhs) / colscale
    return ModeSolution(k=k, a=a, Xp=Xp, roots=roots, anchors=anchors, coeffs=coeffs)


def exact_linear(x: np.ndarray, y: np.ndarray, p: Params, forcing: WindForcing,
                 fields: tuple[str, ...] = ("psi",)) -> dict:
    """Exact steady linear solution on points ``(x, y)``.

    Parameters
    ----------
    x, y:
        Arrays (broadcast together), for example ``grid.X`` and ``grid.Y``.
    p:
        Parameters; ``delta_i`` is ignored (the solution is linear).
    forcing:
        A sine-series wind (any built-in profile). Seasonal modulation is
        ignored (steady solution).
    fields:
        Any of ``"psi", "u", "v", "zeta"``.

    Returns
    -------
    dict mapping field name to array.
    """
    if not forcing.is_sine_series:
        raise ValueError("exact solution needs a sine-series (built-in) wind profile")
    if abs(forcing.Ly - p.Ly) > 1e-12:
        raise ValueError("forcing.Ly and params.Ly differ")
    if p.dm3 > 0 and p.bc_zonal_walls != "free-slip":
        raise ValueError(
            "the separable exact solution needs free-slip southern and northern walls "
            "(set bc_ns='free-slip'); sin(k y) cannot satisfy no slip there"
        )
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    out = {f: np.zeros_like(x) for f in fields}
    for n, a in forcing.modes():
        k = n * math.pi / p.Ly
        m = solve_mode(k, a, p)
        s, c = np.sin(k * y), np.cos(k * y)
        X0 = m.derivative(x, 0)
        if "psi" in out:
            out["psi"] += X0 * s
        if "u" in out:
            out["u"] += -k * X0 * c
        if "v" in out:
            out["v"] += m.derivative(x, 1) * s
        if "zeta" in out:
            out["zeta"] += (m.derivative(x, 2) - k**2 * X0) * s
    return out


# ---------------------------------------------------------------------------
# Sverdrup interior
# ---------------------------------------------------------------------------
def sverdrup_interior(x: np.ndarray, y: np.ndarray, p: Params, forcing: WindForcing) -> np.ndarray:
    """Sverdrup streamfunction ``psi_I = -(1/b) int_x^Lx F dx'``.

    It satisfies ``b psi_x = F`` and ``psi = 0`` on the *eastern* wall; it
    cannot also satisfy the western wall condition, which is why a western
    boundary layer is needed.
    """
    if p.beta_factor == 0:
        raise ValueError("no Sverdrup balance on an f plane (beta_factor = 0)")
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    if forcing.is_sine_series:
        return -(p.Lx - x) * forcing.curl_steady(x, y) / p.beta_factor
    # general x-dependent curl: Gauss-Legendre quadrature from x to the east wall
    nodes, weights = np.polynomial.legendre.leggauss(64)
    half = 0.5 * (p.Lx - x)
    xq = x[..., None] + half[..., None] * (nodes + 1.0)
    F = forcing.curl_steady(xq, np.broadcast_to(y[..., None], xq.shape))
    return -(half * (F * weights).sum(-1)) / p.beta_factor


def sverdrup_transport(y: float, p: Params, forcing: WindForcing) -> float:
    """Sverdrup transport at latitude ``y``: ``psi_I`` at the western edge of the interior."""
    return float(sverdrup_interior(np.array(0.0), np.array(y), p, forcing))


# ---------------------------------------------------------------------------
# boundary-layer theory
# ---------------------------------------------------------------------------
def boundary_layer_profile(x: np.ndarray, p: Params) -> tuple[np.ndarray, np.ndarray]:
    """Classical western-boundary-layer profile, normalized by the interior value.

    Solves ``dM^3 X''' - dS X' - b X = 0`` for the cross-stream structure:
    ``psi_BL(x) = psi_I(0, y) * G(x)`` with ``G(0) = 0``, ``G -> 1`` eastward and,
    when ``delta_m > 0``, ``G''(0) = 0`` (free slip) or ``G'(0) = 0`` (no slip).

    Special cases:

    * Stommel: ``G = 1 - exp(-b x / delta_S)``.
    * Munk, no slip: ``G = 1 - exp(-xi/2) [cos(sqrt3 xi/2) + sin(sqrt3 xi/2)/sqrt3]``, ``xi = x/delta_M``.
    * Munk, free slip: ``G = 1 - exp(-xi/2) [cos(sqrt3 xi/2) - sin(sqrt3 xi/2)/sqrt3]``.

    Returns ``(G, dG/dx)``; ``dG/dx`` is the normalized northward velocity.
    """
    if p.beta_factor <= 0:
        raise ValueError("boundary-layer theory needs beta_factor > 0")
    x = np.asarray(x, float)
    roots = np.array(boundary_layer_roots(p))
    west = roots[roots.real < 0]
    nconds = 2 if p.dm3 > 0 else 1
    if len(west) != nconds:
        raise ValueError("unexpected root structure for the boundary-layer problem")
    rows = [west**0]
    rhs = [-1.0]
    if p.dm3 > 0:
        d = 2 if p.bc == "free-slip" else 1
        rows.append(west**d)
        rhs.append(0.0)
    c = np.linalg.solve(np.array(rows, complex), np.array(rhs, complex))
    e = np.exp(np.multiply.outer(x, west))
    G = 1.0 + (e * c).sum(-1).real
    dG = (e * c * west).sum(-1).real
    return G, dG


def measure_widths(x: np.ndarray, v: np.ndarray) -> dict:
    """Width metrics of a western boundary current from a velocity section.

    ``x`` increases eastward from the western wall; ``v`` is the along-wall
    velocity. The current is taken to be the extremum of ``v`` in the western
    half of the section.

    Returns a dict with

    * ``x_peak``: position of the speed maximum;
    * ``v_peak``: peak velocity;
    * ``e_fold``: distance from the wall to where ``|v|`` first falls below
      ``|v_peak|/e`` east of the peak (linear interpolation);
    * ``zero_crossing``: first sign change of ``v`` east of the peak, i.e.
      the inner edge of the offshore counter-current (``nan`` if none in the
      western half).
    """
    x = np.asarray(x, float)
    v = np.asarray(v, float)
    half = x <= x[0] + 0.5 * (x[-1] - x[0])
    xi, vi = x[half], v[half]
    ip = int(np.argmax(np.abs(vi)))
    vp = vi[ip]
    out = {"x_peak": float(xi[ip]), "v_peak": float(vp), "e_fold": math.nan,
           "zero_crossing": math.nan}
    target = abs(vp) / math.e
    s = np.sign(vp)
    for i in range(ip, len(vi) - 1):
        if s * vi[i + 1] < target <= s * vi[i]:
            f = (s * vi[i] - target) / (s * vi[i] - s * vi[i + 1])
            out["e_fold"] = float(xi[i] + f * (xi[i + 1] - xi[i]) - x[0])
            break
    for i in range(ip, len(vi) - 1):
        if s * vi[i] > 0 >= s * vi[i + 1]:
            f = vi[i] / (vi[i] - vi[i + 1])
            out["zero_crossing"] = float(xi[i] + f * (xi[i + 1] - xi[i]) - x[0])
            break
    return out


def theory_widths(p: Params) -> dict:
    """Boundary-layer width metrics predicted by :func:`boundary_layer_profile`.

    Closed forms: Stommel ``e_fold = delta_S / b``; Munk no-slip
    ``zero_crossing = 2 pi delta_M / sqrt(3)``; Munk free-slip
    ``zero_crossing = 4 pi delta_M / (3 sqrt(3))``.
    """
    roots = np.array(boundary_layer_roots(p))
    west = roots[roots.real < 0]
    scale = 1.0 / float(np.min(np.abs(west.real)))
    xs = np.linspace(0.0, 40.0 * scale, 40001)
    _, dG = boundary_layer_profile(xs, p)
    w = measure_widths(np.concatenate([xs, xs[-1] + xs[1:]]),
                       np.concatenate([dG, np.zeros(len(xs) - 1)]))
    w["decay_scale"] = scale
    w.pop("v_peak", None)
    return w
