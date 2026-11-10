"""Case files: TOML configuration, built-in cases and the case runner.

A case is a TOML file with the sections ``[case]``, ``[grid]``,
``[physics]``, ``[wind]`` and, depending on ``case.kind``, ``[run]``,
``[floats]``, ``[sweep]``, ``[animation]`` and ``[scales]``. Every key has
a default, so a case file only needs the keys it changes. See the
documentation page *Case files* for the full schema.

``kind`` selects what is computed:

``steady``
    The steady linear solution (sparse direct solve), compared with the
    exact separable solution when one exists, with floats circulating
    through the steady gyre in the animation.
``transient``
    Time integration from rest of the full nonlinear equation.
``sweep``
    A sequence of steady solutions while one parameter changes; the
    animation shows the boundary current forming or narrowing.
"""

from __future__ import annotations

import copy
import json
import math
import sys
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any, Optional

import numpy as np

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

from . import diagnostics as diag
from .analytic import exact_linear
from .floats import Floats
from .forcing import WindForcing
from .grid import Grid
from .params import DimensionalScales, Params
from .steady import solve_steady
from .transient import RunConfig, TransientSolver

DEFAULTS: dict[str, dict[str, Any]] = {
    "case": {"name": "custom", "title": "Custom case", "kind": "steady", "description": ""},
    "grid": {"nx": 128, "ny": 128, "Lx": 1.0, "Ly": 1.0},
    "physics": {"delta_s": 0.04, "delta_m": 0.0, "delta_i": 0.0, "beta_factor": 1.0,
                "bc": "free-slip", "bc_ns": None, "from_scales": False},
    "wind": {"profile": "single-gyre", "amplitude": 1.0, "sign": 1, "asymmetry": 0.5,
             "seasonal_amplitude": 0.0, "seasonal_period": 50.0},
    "run": {"t_end": 100.0, "n_frames": 120, "t_record": 0.0, "cfl": 0.6, "dt_max": None,
            "diag_per_frame": 2, "init_noise": 0.0, "init_seed": 0},
    "floats": {"n": 500, "seed": 7, "delta_i": None, "laps": 2.5},
    "sweep": {"parameter": "beta_factor", "start": 0.0, "stop": 1.0, "n": 48, "scale": "linear"},
    "animation": {"fps": 18, "dpi": 85, "trail": 10, "width": 12.0, "height": 7.4},
    "scales": {},
}

KINDS = ("steady", "transient", "sweep")
SWEEPABLE = ("delta_s", "delta_m", "beta_factor", "asymmetry", "amplitude")


# ---------------------------------------------------------------------------
@dataclass
class Case:
    """A fully resolved case (all defaults filled in)."""

    spec: dict
    source: str = "<dict>"

    # convenient accessors ----------------------------------------------
    @property
    def name(self) -> str:
        return self.spec["case"]["name"]

    @property
    def title(self) -> str:
        return self.spec["case"]["title"]

    @property
    def kind(self) -> str:
        return self.spec["case"]["kind"]

    @property
    def description(self) -> str:
        return self.spec["case"]["description"].strip()

    def scales(self) -> Optional[DimensionalScales]:
        sc = self.spec.get("scales") or {}
        if not sc:
            return None
        return DimensionalScales(**sc)

    def grid(self) -> Grid:
        g = self.spec["grid"]
        return Grid(int(g["nx"]), int(g["ny"]), float(g["Lx"]), float(g["Ly"]))

    def params(self, **override) -> Params:
        ph = dict(self.spec["physics"])
        g = self.spec["grid"]
        if ph.pop("from_scales", False):
            sc = self.scales()
            if sc is None:
                raise ValueError("physics.from_scales = true needs a [scales] section")
            ph["delta_s"], ph["delta_m"], ph["delta_i"] = sc.delta_s, sc.delta_m, sc.delta_i
        ph.update({k: v for k, v in override.items() if k in ph})
        return Params(delta_s=float(ph["delta_s"]), delta_m=float(ph["delta_m"]),
                      delta_i=float(ph["delta_i"]), beta_factor=float(ph["beta_factor"]),
                      bc=ph["bc"], bc_ns=ph.get("bc_ns"), Lx=float(g["Lx"]), Ly=float(g["Ly"]))

    def forcing(self, **override) -> WindForcing:
        w = dict(self.spec["wind"])
        w.update({k: v for k, v in override.items() if k in w})
        return WindForcing(profile=w["profile"], amplitude=float(w["amplitude"]),
                           sign=float(w["sign"]), asymmetry=float(w["asymmetry"]),
                           seasonal_amplitude=float(w["seasonal_amplitude"]),
                           seasonal_period=float(w["seasonal_period"]),
                           Ly=float(self.spec["grid"]["Ly"]))

    def run_config(self) -> RunConfig:
        r, fl = self.spec["run"], self.spec["floats"]
        return RunConfig(t_end=float(r["t_end"]), n_frames=int(r["n_frames"]),
                         t_record=float(r["t_record"]), cfl=float(r["cfl"]),
                         dt_max=None if r["dt_max"] is None else float(r["dt_max"]),
                         n_floats=int(fl["n"]), float_seed=int(fl["seed"]),
                         float_delta_i=None if fl["delta_i"] is None else float(fl["delta_i"]),
                         diag_per_frame=int(r["diag_per_frame"]),
                         init_noise=float(r["init_noise"]), init_seed=int(r["init_seed"]))

    def to_toml(self) -> str:
        return dump_toml(self.spec)


# ---------------------------------------------------------------------------
def _parse_value(text: str) -> Any:
    t = text.strip()
    low = t.lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("none", "null"):
        return None
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)
    except ValueError:
        return t.strip("\"'")


def _merge(spec: dict, user: dict) -> dict:
    out = copy.deepcopy(spec)
    for sec, vals in user.items():
        if sec not in out:
            raise ValueError(f"unknown section [{sec}]; expected one of {sorted(out)}")
        if not isinstance(vals, dict):
            raise ValueError(f"section [{sec}] must be a table")
        if sec == "scales":
            out[sec] = dict(vals)
            continue
        for k, v in vals.items():
            if k not in out[sec]:
                raise ValueError(f"unknown key {sec}.{k}; expected one of {sorted(out[sec])}")
            out[sec][k] = v
    return out


def apply_overrides(spec: dict, overrides: list[str] | None) -> dict:
    """Apply ``section.key=value`` strings (from ``--set``)."""
    spec = copy.deepcopy(spec)
    for item in overrides or []:
        if "=" not in item or "." not in item.split("=", 1)[0]:
            raise ValueError(f"override {item!r} must look like section.key=value")
        key, val = item.split("=", 1)
        sec, k = key.strip().split(".", 1)
        if sec not in spec:
            raise ValueError(f"unknown section {sec!r}")
        if sec != "scales" and k not in spec[sec]:
            raise ValueError(f"unknown key {sec}.{k}")
        spec[sec][k] = _parse_value(val)
    return spec


def builtin_names() -> list[str]:
    """Names of the built-in cases."""
    files = resources.files("dewi_lanjar").joinpath("cases")
    return sorted(p.name[:-5] for p in files.iterdir() if p.name.endswith(".toml"))


def load_case(name_or_path: str | Path, overrides: list[str] | None = None) -> Case:
    """Load a built-in case by name or a TOML file by path."""
    p = Path(name_or_path)
    if p.suffix == ".toml" and p.exists():
        text, source = p.read_text(), str(p)
    else:
        name = str(name_or_path)
        files = resources.files("dewi_lanjar").joinpath("cases")
        f = files.joinpath(f"{name}.toml")
        if not f.is_file():
            raise ValueError(f"no built-in case {name!r} and no such file; "
                             f"built-in cases: {', '.join(builtin_names())}")
        text, source = f.read_text(), f"builtin:{name}"
    user = tomllib.loads(text)
    spec = apply_overrides(_merge(DEFAULTS, user), overrides)
    validate_spec(spec)
    return Case(spec=spec, source=source)


def case_from_dict(user: dict, overrides: list[str] | None = None) -> Case:
    spec = apply_overrides(_merge(DEFAULTS, user), overrides)
    validate_spec(spec)
    return Case(spec=spec)


def validate_spec(spec: dict) -> None:
    if spec["case"]["kind"] not in KINDS:
        raise ValueError(f"case.kind must be one of {KINDS}")
    if spec["case"]["kind"] == "sweep" and spec["sweep"]["parameter"] not in SWEEPABLE:
        raise ValueError(f"sweep.parameter must be one of {SWEEPABLE}")
    c = Case(spec)
    c.grid(), c.params(), c.forcing()
    if spec["case"]["kind"] == "transient":
        c.run_config()


def dump_toml(spec: dict) -> str:
    """Minimal TOML writer for case specs (scalars and strings only)."""
    def fmt(v):
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, (int, float)):
            return repr(v)
        s = str(v)
        if "\n" in s:
            return '"""\n' + s.strip("\n") + '\n"""'
        return json.dumps(s)
    lines = []
    for sec, vals in spec.items():
        if not vals:
            continue
        lines.append(f"[{sec}]")
        for k, v in vals.items():
            if v is None:
                continue
            lines.append(f"{k} = {fmt(v)}")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else str(f)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def exact_available(p: Params, forcing: WindForcing) -> bool:
    if not forcing.is_sine_series or not p.has_friction:
        return False
    if p.dm3 > 0 and p.bc_zonal_walls != "free-slip":
        return False
    return True


@dataclass
class CaseResult:
    case: Case
    outdir: Path
    files: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)
    data: dict = field(default_factory=dict)


def run_case(case: Case, outdir: str | Path = "outputs", gif: bool = True, png: bool = True,
             verbose: bool = True) -> CaseResult:
    """Run a case and write its outputs to ``outdir/<case name>/``.

    Outputs: ``summary.txt`` (diagnostics and case description),
    ``diagnostics.json``, ``fields.npz`` (final fields, coordinates and time
    series), ``case.toml`` (fully resolved configuration), a static figure
    ``<name>.png`` and the animation ``<name>.gif``.
    """
    from . import plotting  # lazy: keeps the numerics importable without matplotlib

    out = Path(outdir) / case.name
    out.mkdir(parents=True, exist_ok=True)
    res = CaseResult(case=case, outdir=out)
    log = print if verbose else (lambda *a, **k: None)
    g, p, f = case.grid(), case.params(), case.forcing()
    log(f"[dewi-lanjar] {case.title}")
    log(f"  kind={case.kind}  grid={g.nx}x{g.ny}  {p.regime()['label']}")
    scales = case.scales()

    if case.kind in ("steady", "sweep"):
        p = p.with_(delta_i=0.0)          # the steady solver is linear

    if case.kind == "steady":
        sol = solve_steady(g, p, f)
        d = diag.compute(sol.psi, sol.zeta, g, p, f)
        d["solver_residual"] = sol.residual
        exact = None
        if exact_available(p, f):
            exact = exact_linear(g.X, g.Y, p, f, ("psi", "v"))
            d["exact_max_rel_error"] = float(
                np.max(np.abs(sol.psi - exact["psi"])) / np.max(np.abs(exact["psi"])))
        fl = case.spec["floats"]
        frames_pos = None
        if int(fl["n"]) > 0:
            from .operators import velocity
            u, v = velocity(sol.psi, g.dx, g.dy)
            vmax = float(np.max(np.hypot(u, v)))
            t_total = float(fl["laps"]) * 2 * (g.Lx + g.Ly) / max(vmax, 1e-12)
            n_frames = int(case.spec["run"]["n_frames"])
            frames_pos = Floats(g, int(fl["n"]), seed=int(fl["seed"])).advect_steady(
                sol.psi, t_total, n_frames)
        res.data = {"psi": sol.psi, "zeta": sol.zeta, "exact": exact, "floats": frames_pos}
        res.diagnostics = d
        if png:
            res.files["png"] = plotting.steady_figure(case, sol.psi, sol.zeta, d, exact,
                                                      out / f"{case.name}.png", scales)
        if gif and frames_pos is not None:
            res.files["gif"] = plotting.animate_steady(case, sol.psi, sol.zeta, d, exact,
                                                       frames_pos, out / f"{case.name}.gif")
        np.savez_compressed(out / "fields.npz", x=g.x, y=g.y, psi=sol.psi, zeta=sol.zeta,
                            **({"psi_exact": exact["psi"]} if exact else {}))

    elif case.kind == "transient":
        cfg = case.run_config()
        solver = TransientSolver(g, p, f)
        last = [-1]

        def progress(t, T):
            pct = int(100 * t / T)
            if verbose and pct // 10 != last[0]:
                last[0] = pct // 10
                log(f"  t = {t:8.2f} / {T:g}  ({pct:3d}%)")

        r = solver.run(cfg, progress=progress)
        d = r.final
        d["stats"] = r.stats
        res.data = {"result": r}
        res.diagnostics = d
        if png:
            res.files["png"] = plotting.transient_figure(case, r, out / f"{case.name}.png", scales)
        if gif:
            res.files["gif"] = plotting.animate_transient(case, r, out / f"{case.name}.gif")
        np.savez_compressed(out / "fields.npz", x=g.x, y=g.y, times=r.times,
                            psi_final=r.psi[-1], zeta_final=r.zeta[-1],
                            **{f"series_{k}": v for k, v in r.series.items()})

    else:  # sweep
        sw = case.spec["sweep"]
        if sw["scale"] == "log":
            values = np.geomspace(float(sw["start"]), float(sw["stop"]), int(sw["n"]))
        else:
            values = np.linspace(float(sw["start"]), float(sw["stop"]), int(sw["n"]))
        par = sw["parameter"]
        frames, dlist = [], []
        for val in values:
            pp = (case.params(**{par: float(val)}) if par in case.spec["physics"] else p).with_(
                delta_i=0.0)
            ff = case.forcing(**{par: float(val)}) if par in case.spec["wind"] else f
            sol = solve_steady(g, pp, ff)
            dd = diag.compute(sol.psi, sol.zeta, g, pp, ff, full=True)
            dd["value"] = float(val)
            frames.append(sol.psi)
            dlist.append(dd)
        d = {"parameter": par, "values": values.tolist(),
             "intensification_index": [x["intensification_index"] for x in dlist],
             "ke_centroid_x": [x["ke_centroid_x"] for x in dlist],
             "ke_zonal_entropy": [x["ke_zonal_entropy"] for x in dlist],
             "wbc_e_fold": [x["wbc_e_fold"] for x in dlist],
             "theory_e_fold": [x["theory_e_fold"] for x in dlist],
             "final": dlist[-1]}
        res.data = {"values": values, "psi": np.array(frames), "diags": dlist}
        res.diagnostics = d
        if png:
            res.files["png"] = plotting.sweep_figure(case, values, dlist, out / f"{case.name}.png")
        if gif:
            res.files["gif"] = plotting.animate_sweep(case, values, np.array(frames), dlist,
                                                      out / f"{case.name}.gif")
        np.savez_compressed(out / "fields.npz", x=g.x, y=g.y, values=values, psi=np.array(frames))

    # text outputs ----------------------------------------------------
    summary = [f"dewi-lanjar case: {case.title}", f"source: {case.source}", ""]
    if case.description:
        summary += [case.description, ""]
    final = res.diagnostics["final"] if case.kind == "sweep" else res.diagnostics
    summary.append(diag.describe(final, scales))
    if case.kind == "steady":
        summary.append("")
        summary.append(f"Sparse solver relative residual: {final['solver_residual']:.2e}")
        if "exact_max_rel_error" in final:
            summary.append("Max relative error against the exact separable solution: "
                           f"{final['exact_max_rel_error']:.3e}  (second order in h)")
    if case.kind == "transient":
        st = final["stats"]
        summary.append("")
        summary.append(f"Time steps: {st['steps']}, dt in [{st['dt_min']:.3g}, {st['dt_max']:.3g}], "
                       f"wall time {st['wall_seconds']:.1f} s")
    if scales is not None:
        summary += ["", scales.explain()]
    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    (out / "diagnostics.json").write_text(json.dumps(_jsonable(res.diagnostics), indent=2))
    (out / "case.toml").write_text(case.to_toml())
    res.files.update({"summary": out / "summary.txt", "json": out / "diagnostics.json",
                      "npz": out / "fields.npz", "toml": out / "case.toml"})
    log(f"  wrote {out}")
    return res
