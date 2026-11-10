"""Static figures and GIF animations.

All rendering uses the non-interactive Agg backend (no GUI). Animations
are written frame by frame and encoded with a single global palette and no
dithering, which keeps GIFs small and free of flicker.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from PIL import Image  # noqa: E402

from .analytic import boundary_layer_profile  # noqa: E402

# ---------------------------------------------------------------------------
# palette
BG = "#0a0f1e"
PANEL = "#0d1428"
FG = "#e6ecff"
MUTED = "#8a96b8"
ACCENT = "#ffd166"
CYAN = "#5ee6ff"
PINK = "#ff5c8a"

PSI_CMAP = LinearSegmentedColormap.from_list(
    "dl_psi", ["#9ff7ff", "#2aa3e0", "#123f7a", "#0a0f1e", "#5e1640", "#e2457a", "#ffd88a"], N=512)
ZETA_CMAP = LinearSegmentedColormap.from_list(
    "dl_zeta", ["#c9f5ff", "#3fb6ff", "#1d3f8f", "#0a0f1e", "#7a2a12", "#ff8a3d", "#fff0b0"], N=512)


def _style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
        "axes.edgecolor": "#2a3557", "axes.labelcolor": FG, "text.color": FG,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": FG,
        "font.size": 9, "axes.titlesize": 10, "font.family": "DejaVu Sans",
        "axes.titlepad": 6, "axes.labelsize": 9,
        "axes.grid": False, "legend.frameon": False, "legend.fontsize": 8,
    })


def _asinh_norm(z: np.ndarray, frac: float = 0.03) -> tuple[np.ndarray, float]:
    """Compress vorticity with ``asinh(z/z0)`` so boundary layer and interior both show."""
    zmax = float(np.max(np.abs(z))) or 1.0
    z0 = frac * zmax
    return np.arcsinh(z / z0), float(np.arcsinh(zmax / z0))


def _img(ax, field, grid, cmap, vlim):
    return ax.imshow(field.T, origin="lower", extent=(0, grid.Lx, 0, grid.Ly), cmap=cmap,
                     vmin=-vlim, vmax=vlim, interpolation="bilinear", aspect="equal")


def _decorate_map(ax, title):
    ax.set_title(title, loc="left")
    ax.set_xlabel("x / L   (west → east)")
    ax.set_ylabel("y / L   (south → north)")
    ax.tick_params(length=2, labelsize=8)


def _psi_levels(vlim, n=14):
    lv = np.linspace(-vlim, vlim, n + 1)
    return lv[np.abs(lv) > 1e-12 * vlim]


# ---------------------------------------------------------------------------
# GIF encoding
class GifWriter:
    """Collect RGB frames and write an optimized GIF with one global palette."""

    def __init__(self, path: Path, fps: int = 18):
        self.path = Path(path)
        self.fps = fps
        self.frames: list[np.ndarray] = []

    def add(self, fig) -> None:
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        self.frames.append(buf)

    def close(self) -> Path:
        if not self.frames:
            raise RuntimeError("no frames")
        idx = np.linspace(0, len(self.frames) - 1, min(6, len(self.frames))).astype(int)
        sample = Image.fromarray(np.concatenate([self.frames[i] for i in idx], axis=0))
        pal = sample.quantize(colors=255, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
        imgs = [Image.fromarray(f).quantize(palette=pal, dither=Image.Dither.NONE) for f in self.frames]
        dur = int(round(1000 / self.fps))
        imgs[0].save(self.path, save_all=True, append_images=imgs[1:], duration=dur, loop=0,
                     optimize=True, disposal=1)
        self.frames.clear()
        return self.path


def _floats_rgba(history: np.ndarray, f: int, trail: int, color=(0.92, 0.97, 1.0)):
    pts, cols, sizes = [], [], []
    for k in range(trail, -1, -1):
        idx = f - k
        if idx < 0:
            continue
        a = (1 - k / (trail + 1)) ** 1.6
        pts.append(history[idx])
        c = np.empty((history.shape[1], 4))
        c[:, :3] = color
        c[:, 3] = 0.7 * a
        cols.append(c)
        sizes.append(np.full(history.shape[1], 0.3 + 3.4 * a**2))
    return np.concatenate(pts), np.concatenate(cols), np.concatenate(sizes)


# ---------------------------------------------------------------------------
# layout helpers: every label, legend and colorbar lives outside the data
LEGEND_KW = dict(frameon=False, fontsize=8.5, handlelength=1.8, borderaxespad=0.0,
                 labelcolor=FG)
PRETTY = {"beta_factor": "β factor  b", "delta_s": "Stommel width  δ_S",
          "delta_m": "Munk width  δ_M", "asymmetry": "wind asymmetry", "amplitude": "wind amplitude"}

SHORT = {"beta_factor": "b", "delta_s": "δ_S", "delta_m": "δ_M", "asymmetry": "asymmetry",
         "amplitude": "amplitude"}


def _legend_right(ax, **kw):
    """Legend in a column to the right of the axes, outside the plotting area."""
    return ax.legend(loc="center left", bbox_to_anchor=(1.03, 0.5), **{**LEGEND_KW, **kw})


def _colorbar(fig, ax, mappable, label):
    """Colorbar attached to the displayed map (respects equal aspect), label outside."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.08)
    cb = fig.colorbar(mappable, cax=cax)
    cb.outline.set_edgecolor("#2a3557")
    cb.ax.tick_params(labelsize=7.5, length=2, colors=MUTED)
    cb.set_label(label, color=FG, fontsize=9)
    return cb


def _header(fig, title, subtitle="", live=False):
    """Title, a static subtitle and (optionally) a live readout line, above all panels."""
    fig.text(0.05, 0.962, title, fontsize=13.5, color=FG, va="center")
    if subtitle:
        fig.text(0.05, 0.928, subtitle, fontsize=9, color=MUTED, va="center")
    fig.text(0.985, 0.012, "dewi-lanjar 0.0.1", fontsize=7.5, color="#4a5578", ha="right")
    if live:
        return fig.text(0.05, 0.895, "", fontsize=9, color=ACCENT, va="center",
                        family="DejaVu Sans Mono")
    return None


def _map_axes(fig):
    """Two map panels in the top row (the colorbars are appended to their right)."""
    ax1 = fig.add_axes([0.05, 0.40, 0.36, 0.43])
    ax2 = fig.add_axes([0.53, 0.40, 0.36, 0.43])
    for ax in (ax1, ax2):
        ax.set_anchor("W")
    return ax1, ax2


def _subtitle(case, p):
    reg = p.regime()
    parts = [reg["label"],
             f"δ_S = {p.delta_s:.3g}", f"δ_M = {p.delta_m:.3g}", f"δ_I = {p.delta_i:.3g}",
             f"b = {p.beta_factor:g}", f"walls: {p.bc} / {p.bc_zonal_walls}"]
    return "   ·   ".join(parts)


def _summary_line(d):
    parts = []
    if math.isfinite(d.get("intensification_index", math.nan)):
        parts.append(f"intensification {d['intensification_index']:.1f}")
    if math.isfinite(d.get("wbc_e_fold", math.nan)):
        s = f"e-fold width {d['wbc_e_fold']:.3f}"
        if math.isfinite(d.get("theory_e_fold", math.nan)):
            s += f" (theory {d['theory_e_fold']:.3f})"
        parts.append(s)
    if math.isfinite(d.get("transport_ratio", math.nan)):
        parts.append(f"T / T_Sverdrup {d['transport_ratio']:.2f}")
    parts.append(f"KE centroid x/L {d['ke_centroid_x']:.2f}")
    return "   ·   ".join(parts)


# ---------------------------------------------------------------------------
# panels
def _section_panel(ax, grid, p, psi, d, exact=None, forcing=None):
    from .operators import velocity
    j = int(np.argmin(np.abs(grid.y - d["section_y"])))
    _, v = velocity(psi, grid.dx, grid.dy)
    ax.axhline(0, color="#2a3557", lw=0.8)
    ax.plot(grid.x, v[:, j], color=CYAN, lw=1.8, label="numerical")
    if exact is not None:
        ax.plot(grid.x, exact["v"][:, j], color=ACCENT, lw=1.0, ls="--", label="exact separable")
    if p.beta_factor > 0 and p.has_friction and forcing is not None and forcing.is_sine_series:
        # uniformly valid composite psi = psi_I(x, y) G(x); psi_I = -(Lx - x) F(y) / b
        G, dG = boundary_layer_profile(grid.x, p)
        c = -float(forcing.curl_steady(0.0, d["section_y"])) / p.beta_factor
        ax.plot(grid.x, -c * G + (p.Lx - grid.x) * c * dG, color=PINK, lw=1.0, ls=":",
                label="Sverdrup + layer theory")
    if math.isfinite(d.get("wbc_e_fold", math.nan)):
        ax.axvline(d["wbc_e_fold"], color=MUTED, lw=0.8, ls="--", label="e-folding width")
    ax.set_xlim(0, grid.Lx)
    ax.set_xlabel("x / L")
    ax.set_ylabel("v (northward)")
    ax.set_title(f"northward velocity along y = {d['section_y']:.2f}", loc="left")
    _legend_right(ax)


def _budget_panel(ax, d):
    vb = d.get("vorticity_budget")
    if not vb:
        ax.axis("off")
        return
    from matplotlib.patches import Patch
    terms = ["wind", "beta", "drag", "viscosity", "advection"]
    cols = [ACCENT, CYAN, "#9b8cff", PINK, "#7dffb0"]
    xw = np.arange(len(terms))
    w = 0.38
    ax.bar(xw - w / 2, [vb["west"][t] for t in terms], w, color=cols, alpha=0.95)
    ax.bar(xw + w / 2, [vb["interior"][t] for t in terms], w, color=cols, alpha=0.4,
           hatch="///", edgecolor=BG)
    ax.axhline(0, color="#2a3557", lw=0.8)
    ax.set_xticks(xw, terms, rotation=20, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", labelsize=8)
    ax.set_ylabel("integrated tendency")
    ax.set_title("vorticity budget of the gyre", loc="left")
    handles = [Patch(facecolor=MUTED, alpha=0.95, label="western strip"),
               Patch(facecolor=MUTED, alpha=0.4, hatch="///", edgecolor=BG, label="interior")]
    _legend_right(ax, handles=handles)


def _bottom_pair(fig):
    """Two bottom panels, each with a legend column to its right."""
    return fig.add_axes([0.06, 0.09, 0.27, 0.21]), fig.add_axes([0.57, 0.09, 0.27, 0.21])


def _bottom_wide(fig):
    """One wide bottom panel with a legend column to its right."""
    return fig.add_axes([0.06, 0.08, 0.66, 0.22])


def _series_panel(ax, s, t_end):
    ke_n = s["KE"] / (np.max(s["KE"]) or 1.0)
    tr = np.asarray(s["transport_ratio"], float)
    ax.plot(s["t"], ke_n, color=CYAN, lw=1.5, label="KE / max KE")
    if np.any(np.isfinite(tr)):
        ax.plot(s["t"], tr, color=ACCENT, lw=1.2, label="WBC transport / Sverdrup")
    ax.plot(s["t"], s["ke_centroid_x"], color=PINK, lw=1.1, label="KE zonal centroid x/L")
    ax.plot(s["t"], s["ke_zonal_entropy"], color="#9b8cff", lw=1.1, label="zonal KE entropy")
    ax.axhline(1, color="#2a3557", lw=0.6, ls="--")
    ax.set_xlim(0, t_end)
    finite = np.concatenate([ke_n, tr[np.isfinite(tr)], s["ke_centroid_x"], s["ke_zonal_entropy"]])
    ax.set_ylim(min(0.0, float(np.nanmin(finite))) - 0.05,
                min(4.0, max(1.2, float(np.nanpercentile(finite, 99.5)) * 1.1)))
    ax.set_xlabel("t  [1/(βL)]")
    ax.set_title("gyre diagnostics through time", loc="left")
    _legend_right(ax)


def _energy_panel(ax, s):
    ax.plot(s["t"], s["energy_wind"], color=ACCENT, lw=1.2, label="wind work W")
    ax.plot(s["t"], s["energy_drag"], color="#9b8cff", lw=1.2, label="drag")
    ax.plot(s["t"], s["energy_viscous"], color=PINK, lw=1.2, label="viscous")
    ax.plot(s["t"], s["dKE_dt"], color=FG, lw=0.8, ls="--", label="dKE/dt")
    ax.axhline(0, color="#2a3557", lw=0.8)
    ax.set_xlabel("t  [1/(βL)]")
    ax.set_title("energy budget", loc="left")
    _legend_right(ax)


# ---------------------------------------------------------------------------
# static figures
def _static_fig():
    _style()
    return plt.figure(figsize=(12, 7.4))


def steady_figure(case, psi, zeta, d, exact, path, scales=None):
    g, p = case.grid(), case.params().with_(delta_i=0.0)     # the steady solver is linear
    fig = _static_fig()
    _header(fig, case.title, _subtitle(case, p))
    fig.text(0.05, 0.895, _summary_line(d), fontsize=9, color=ACCENT, va="center")
    ax1, ax2 = _map_axes(fig)
    vlim = float(np.max(np.abs(psi))) or 1.0
    im = _img(ax1, psi, g, PSI_CMAP, vlim)
    ax1.contour(g.x, g.y, psi.T, levels=_psi_levels(vlim), colors="w", linewidths=0.5, alpha=0.4)
    _decorate_map(ax1, "streamfunction ψ")
    _colorbar(fig, ax1, im, "ψ")
    zz, zl = _asinh_norm(zeta)
    im2 = _img(ax2, zz, g, ZETA_CMAP, zl)
    _decorate_map(ax2, "relative vorticity ζ")
    _colorbar(fig, ax2, im2, "asinh(ζ / ζ₀)")
    b1, b2 = _bottom_pair(fig)
    _section_panel(b1, g, p, psi, d, exact, case.forcing())
    _budget_panel(b2, d)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return Path(path)


def transient_figure(case, r, path, scales=None):
    g = r.grid
    psi, zeta = r.psi[-1].astype(float), r.zeta[-1].astype(float)
    fig = _static_fig()
    _header(fig, case.title, _subtitle(case, r.params))
    fig.text(0.05, 0.895, f"t = {r.times[-1]:.1f}   ·   " + _summary_line(r.final),
             fontsize=9, color=ACCENT, va="center")
    ax1, ax2 = _map_axes(fig)
    vlim = float(np.max(np.abs(psi))) or 1.0
    im = _img(ax1, psi, g, PSI_CMAP, vlim)
    ax1.contour(g.x, g.y, psi.T, levels=_psi_levels(vlim), colors="w", linewidths=0.5, alpha=0.4)
    _decorate_map(ax1, "streamfunction ψ (final state)")
    _colorbar(fig, ax1, im, "ψ")
    zz, zl = _asinh_norm(zeta)
    im2 = _img(ax2, zz, g, ZETA_CMAP, zl)
    _decorate_map(ax2, "relative vorticity ζ (final state)")
    _colorbar(fig, ax2, im2, "asinh(ζ / ζ₀)")
    b1, b2 = _bottom_pair(fig)
    _series_panel(b1, r.series, r.times[-1])
    _energy_panel(b2, r.series)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return Path(path)


def sweep_figure(case, values, dlist, path):
    _style()
    par = case.spec["sweep"]["parameter"]
    fig = plt.figure(figsize=(12, 4.2))
    fig.text(0.05, 0.94, case.title, fontsize=13, color=FG, va="center")
    axs = [fig.add_axes([0.06 + k * 0.32, 0.30, 0.24, 0.52]) for k in range(3)]
    axs[0].plot(values, [d["intensification_index"] for d in dlist], color=ACCENT)
    axs[0].set_title("intensification index", loc="left")
    axs[1].plot(values, [d["ke_centroid_x"] for d in dlist], color=CYAN, label="KE centroid x/L")
    axs[1].plot(values, [d["ke_zonal_entropy"] for d in dlist], color=PINK, label="zonal KE entropy")
    axs[1].set_title("energy concentration", loc="left")
    axs[2].plot(values, [d["wbc_e_fold"] for d in dlist], color=CYAN, label="measured")
    axs[2].plot(values, [d["theory_e_fold"] for d in dlist], color=ACCENT, ls="--", label="theory")
    axs[2].set_title("e-folding width", loc="left")
    for a in axs:
        a.set_xlabel(PRETTY.get(par, par))
        if case.spec["sweep"]["scale"] == "log":
            a.set_xscale("log")
    for a in axs[1:]:
        a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, **LEGEND_KW)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return Path(path)


# ---------------------------------------------------------------------------
# animations
def _anim_figure(case):
    a = case.spec["animation"]
    _style()
    return plt.figure(figsize=(float(a["width"]), float(a["height"])), dpi=int(a["dpi"]))


def _set_floats(sc, history, f, trail):
    pts, cols, sz = _floats_rgba(history, f, trail)
    sc.set_offsets(pts)
    sc.set_facecolors(cols)
    sc.set_sizes(sz)


def animate_steady(case, psi, zeta, d, exact, history, path):
    g, p = case.grid(), case.params().with_(delta_i=0.0)     # the steady solver is linear
    a = case.spec["animation"]
    trail = int(a["trail"])
    fig = _anim_figure(case)
    _header(fig, case.title, _subtitle(case, p))
    fig.text(0.05, 0.895, _summary_line(d), fontsize=9, color=ACCENT, va="center")
    ax1, ax2 = _map_axes(fig)
    vlim = float(np.max(np.abs(psi))) or 1.0
    im = _img(ax1, psi, g, PSI_CMAP, vlim)
    ax1.contour(g.x, g.y, psi.T, levels=_psi_levels(vlim), colors="w", linewidths=0.4, alpha=0.35)
    _decorate_map(ax1, "streamfunction ψ with Lagrangian floats")
    _colorbar(fig, ax1, im, "ψ")
    zz, zl = _asinh_norm(zeta)
    im2 = _img(ax2, zz, g, ZETA_CMAP, zl)
    _decorate_map(ax2, "relative vorticity ζ with floats")
    _colorbar(fig, ax2, im2, "asinh(ζ / ζ₀)")
    b1, b2 = _bottom_pair(fig)
    _section_panel(b1, g, p, psi, d, exact, case.forcing())
    _budget_panel(b2, d)
    pts, cols, sz = _floats_rgba(history, 0, trail)
    scs = [ax.scatter(pts[:, 0], pts[:, 1], s=sz, c=cols, linewidths=0) for ax in (ax1, ax2)]
    for ax in (ax1, ax2):
        ax.set_xlim(0, g.Lx)
        ax.set_ylim(0, g.Ly)
    w = GifWriter(path, int(a["fps"]))
    for f in range(len(history)):
        for sc in scs:
            _set_floats(sc, history, f, trail)
        w.add(fig)
    plt.close(fig)
    return w.close()


def animate_transient(case, r, path):
    g = r.grid
    a = case.spec["animation"]
    trail = int(a["trail"])
    fig = _anim_figure(case)
    info = _header(fig, case.title, _subtitle(case, r.params), live=True)
    ax1, ax2 = _map_axes(fig)
    psi_all = r.psi.astype(float)
    vlim = float(np.max(np.abs(psi_all))) or 1.0
    levels = _psi_levels(vlim)
    z0 = 0.03 * (float(np.max(np.abs(r.zeta[-1].astype(float)))) or 1.0)
    zl = float(np.arcsinh(np.max(np.abs(r.zeta.astype(float))) / z0))
    im1 = _img(ax1, psi_all[0], g, PSI_CMAP, vlim)
    im2 = _img(ax2, np.arcsinh(r.zeta[0].astype(float) / z0), g, ZETA_CMAP, zl)
    _decorate_map(ax1, "streamfunction ψ" + (" with Lagrangian floats" if r.floats is not None else ""))
    _decorate_map(ax2, "relative vorticity ζ")
    _colorbar(fig, ax1, im1, "ψ")
    _colorbar(fig, ax2, im2, "asinh(ζ / ζ₀)")
    for ax in (ax1, ax2):
        ax.set_xlim(0, g.Lx)
        ax.set_ylim(0, g.Ly)
    ax3 = _bottom_wide(fig)
    s = r.series
    _series_panel(ax3, s, r.times[-1])
    cursor = ax3.axvline(0, color=FG, lw=1.0, alpha=0.8)
    sc1 = None
    if r.floats is not None:
        pts, cols, sz = _floats_rgba(r.floats, 0, trail)
        sc1 = ax1.scatter(pts[:, 0], pts[:, 1], s=sz, c=cols, linewidths=0)
    contour = [None]
    w = GifWriter(path, int(a["fps"]))
    ts = np.asarray(s["t"])
    for f in range(len(r.times)):
        im1.set_data(psi_all[f].T)
        im2.set_data(np.arcsinh(r.zeta[f].astype(float) / z0).T)
        if contour[0] is not None:
            contour[0].remove()
        contour[0] = ax1.contour(g.x, g.y, psi_all[f].T, levels=levels, colors="w",
                                 linewidths=0.4, alpha=0.35)
        if sc1 is not None:
            _set_floats(sc1, r.floats, f, trail)
        t = r.times[f]
        cursor.set_xdata([t, t])
        k = int(np.argmin(np.abs(ts - t)))
        info.set_text(f"t = {t:<7.1f}  KE = {s['KE'][k]:<9.4g}  "
                      f"intensification = {s['intensification_index'][k]:<6.1f}  "
                      f"ψ ∈ [{s['psi_min'][k]:.2f}, {s['psi_max'][k]:.2f}]")
        w.add(fig)
    plt.close(fig)
    return w.close()


def animate_sweep(case, values, frames, dlist, path):
    from .operators import velocity
    g = case.grid()
    a = case.spec["animation"]
    par = case.spec["sweep"]["parameter"]
    fig = _anim_figure(case)
    p0 = case.params()
    sub = "   ·   ".join([f"sweep of {PRETTY.get(par, par)}",
                          f"from {values[0]:.4g} to {values[-1]:.4g}",
                          f"walls: {p0.bc} / {p0.bc_zonal_walls}"])
    info = _header(fig, case.title, sub, live=True)
    ax1 = fig.add_axes([0.05, 0.40, 0.36, 0.43])
    ax1.set_anchor("W")
    ax2 = fig.add_axes([0.53, 0.40, 0.26, 0.43])
    j = int(np.argmin(np.abs(g.y - dlist[-1]["section_y"])))
    vs = np.array([velocity(f, g.dx, g.dy)[1][:, j] for f in frames])
    vmax = float(np.max(np.abs(vs))) or 1.0
    im = _img(ax1, frames[0] / (np.max(np.abs(frames[0])) or 1), g, PSI_CMAP, 1.0)
    _decorate_map(ax1, "streamfunction ψ / max|ψ|")
    _colorbar(fig, ax1, im, "ψ / max|ψ|")
    line, = ax2.plot(g.x, vs[0], color=CYAN, lw=1.8, label="v at current value")
    ax2.axhline(0, color="#2a3557", lw=0.8)
    efold = ax2.axvline(0, color=ACCENT, lw=0.9, ls="--", label="e-folding width")
    ax2.set_xlim(0, g.Lx)
    ax2.set_ylim(-1.05 * vmax, 1.05 * vmax)
    ax2.set_xlabel("x / L")
    ax2.set_ylabel("v (northward)")
    ax2.set_title(f"northward velocity along y = {g.y[j]:.2f}", loc="left")
    _legend_right(ax2)
    ax3 = _bottom_wide(fig)
    xs = values
    ii = np.array([d["intensification_index"] for d in dlist])
    ax3.plot(xs, ii / np.nanmax(ii), color=ACCENT, lw=1.4, label="intensification (normalized)")
    ax3.plot(xs, [d["ke_centroid_x"] for d in dlist], color=PINK, lw=1.2, label="KE zonal centroid x/L")
    ax3.plot(xs, [d["ke_zonal_entropy"] for d in dlist], color="#9b8cff", lw=1.2,
             label="zonal KE entropy")
    if case.spec["sweep"]["scale"] == "log":
        ax3.set_xscale("log")
    ax3.set_xlim(min(xs), max(xs))
    ax3.set_xlabel(PRETTY.get(par, par))
    ax3.set_title("diagnostics across the sweep", loc="left")
    mk, = ax3.plot([], [], "o", color=FG, ms=6, clip_on=False, zorder=5)
    cur = ax3.axvline(xs[0], color=FG, lw=0.8, alpha=0.6)
    _legend_right(ax3)
    w = GifWriter(path, int(a["fps"]))
    order = list(range(len(values))) + [len(values) - 1] * 8
    for f in order:
        fr = frames[f]
        im.set_data((fr / (np.max(np.abs(fr)) or 1)).T)
        line.set_ydata(vs[f])
        d = dlist[f]
        if math.isfinite(d["wbc_e_fold"]):
            efold.set_xdata([d["wbc_e_fold"]] * 2)
        mk.set_data([xs[f]], [d["ke_centroid_x"]])
        cur.set_xdata([xs[f]] * 2)
        info.set_text(f"{SHORT.get(par, par)} = {xs[f]:<8.4g}  "
                      f"intensification = {d['intensification_index']:<7.2f}  "
                      f"KE centroid x/L = {d['ke_centroid_x']:.3f}")
        w.add(fig)
    plt.close(fig)
    return w.close()
