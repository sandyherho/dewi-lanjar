"""Command-line interface: ``dewi-lanjar <command> ...``.

Commands
--------
list        list the built-in cases
show        print a case (fully resolved TOML) and its description
init        write a case file to edit, starting from a built-in case
run         run a case (built-in name or TOML file) and write figures, GIF, data
run-all     run every built-in case
steady      solve a steady linear gyre from command-line parameters
scales      boundary-layer widths and scales from dimensional ocean parameters
widths      theoretical boundary-layer widths for nondimensional parameters
validate    convergence of the solvers against the exact solution
explain     print a derivation (overview, stommel, munk, ...)
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from ._version import __version__


def _cmd_list(args) -> int:
    from .cases import builtin_names, load_case
    print(f"{'case':<26s}{'kind':<11s}title")
    for n in builtin_names():
        c = load_case(n)
        print(f"{n:<26s}{c.kind:<11s}{c.title}")
    return 0


def _cmd_show(args) -> int:
    from .cases import load_case
    c = load_case(args.case, args.set)
    print(f"# {c.title}\n# source: {c.source}\n")
    if c.description:
        print("\n".join("# " + ln for ln in c.description.splitlines()))
        print()
    print(c.to_toml())
    p = c.params()
    print(f"# regime: {p.regime()['label']}")
    return 0


def _cmd_init(args) -> int:
    from .cases import load_case
    c = load_case(args.template, args.set)
    c.spec["case"]["name"] = Path(args.file).stem
    path = Path(args.file)
    if path.exists() and not args.force:
        print(f"{path} exists; use --force to overwrite", file=sys.stderr)
        return 1
    path.write_text(c.to_toml())
    print(f"wrote {path} (from {args.template}); edit it and run: dewi-lanjar run {path}")
    return 0


def _cmd_run(args) -> int:
    from .cases import load_case, run_case
    c = load_case(args.case, args.set)
    res = run_case(c, args.out, gif=not args.no_gif, png=not args.no_png, verbose=not args.quiet)
    if not args.quiet:
        print()
        print(res.files["summary"].read_text())
    return 0


def _cmd_run_all(args) -> int:
    from .cases import builtin_names, load_case, run_case
    kinds = set(args.kinds.split(",")) if args.kinds else None
    for n in builtin_names():
        c = load_case(n, args.set)
        if kinds and c.kind not in kinds:
            continue
        run_case(c, args.out, gif=not args.no_gif, png=not args.no_png, verbose=not args.quiet)
    return 0


def _params_from_args(args):
    from .forcing import WindForcing
    from .grid import Grid
    from .params import Params
    p = Params(delta_s=args.delta_s, delta_m=args.delta_m, beta_factor=args.beta_factor,
               bc=args.bc, bc_ns=args.bc_ns, Lx=1.0, Ly=args.Ly)
    f = WindForcing(profile=args.wind, sign=args.sign, asymmetry=args.asymmetry, Ly=args.Ly)
    g = Grid(args.nx, args.ny if args.ny else int(round(args.nx * args.Ly)), 1.0, args.Ly)
    return g, p, f


def _cmd_steady(args) -> int:
    import numpy as np

    from . import diagnostics as diag
    from .analytic import exact_linear
    from .cases import exact_available
    from .steady import solve_steady
    g, p, f = _params_from_args(args)
    sol = solve_steady(g, p, f)
    d = diag.compute(sol.psi, sol.zeta, g, p, f)
    print(f"Wind: {f.describe()}")
    print(f"Grid: {g.nx} x {g.ny}, residual {sol.residual:.2e}")
    print(diag.describe(d))
    if exact_available(p, f):
        ex = exact_linear(g.X, g.Y, p, f)["psi"]
        print(f"\nMax relative error vs exact separable solution: "
              f"{np.max(np.abs(sol.psi - ex)) / np.max(np.abs(ex)):.3e}")
    if args.why:
        from .explain import explain
        topic = "fplane" if p.beta_factor == 0 else ("munk" if p.delta_m > 0 else "stommel")
        print()
        print(explain(topic))
    if args.png:
        from .cases import case_from_dict
        from .plotting import steady_figure
        case = case_from_dict({"case": {"name": Path(args.png).stem, "title": "Steady gyre"},
                               "grid": {"nx": g.nx, "ny": g.ny, "Ly": args.Ly}})
        case.spec["physics"].update({"delta_s": p.delta_s, "delta_m": p.delta_m,
                                     "beta_factor": p.beta_factor, "bc": p.bc, "bc_ns": p.bc_ns})
        case.spec["wind"].update({"profile": f.profile, "sign": f.sign, "asymmetry": f.asymmetry})
        exact = exact_linear(g.X, g.Y, p, f, ("psi", "v")) if exact_available(p, f) else None
        steady_figure(case, sol.psi, sol.zeta, d, exact, args.png)
        print(f"\nwrote {args.png}")
    return 0


def _cmd_scales(args) -> int:
    from .params import DimensionalScales, beta_at_latitude
    beta = args.beta if args.beta is not None else beta_at_latitude(args.lat)
    s = DimensionalScales(L=args.L * 1e3, beta=beta, r=args.r, A_h=args.A_h, tau0=args.tau0,
                          rho0=args.rho0, H=args.H)
    if args.beta is None:
        print(f"beta at {args.lat:g} deg latitude = 2 Omega cos(phi)/a = {beta:.4g} 1/(m s)\n")
    print(s.explain())
    p = s.to_params()
    print("\nRegime:", p.regime()["label"])
    return 0


def _cmd_widths(args) -> int:
    from .analytic import theory_widths
    from .params import Params, boundary_layer_roots
    p = Params(delta_s=args.delta_s, delta_m=args.delta_m, beta_factor=args.beta_factor, bc=args.bc)
    roots = boundary_layer_roots(p)
    print("Local boundary-layer equation: dM^3 m^3 - dS m - b = 0")
    print(f"  dS = {p.delta_s:g}, dM = {p.delta_m:g}, b = {p.beta_factor:g}, bc = {p.bc}")
    for r in roots:
        tag = "decays eastward (western layer)" if r.real < 0 else "grows eastward (excluded)"
        print(f"  m = {r.real: .6g} {r.imag:+.6g} i   {tag}")
    w = theory_widths(p)
    print(f"\nDecay scale 1/|Re m|        {w['decay_scale']:.6g}")
    print(f"Peak velocity at x          {w['x_peak']:.6g}")
    print(f"e-folding width             {w['e_fold']:.6g}")
    zc = w["zero_crossing"]
    print(f"Counter-current starts at   {zc:.6g}" if math.isfinite(zc) else
          "Counter-current             none (monotone decay)")
    return 0


def _cmd_validate(args) -> int:
    import numpy as np

    from .analytic import exact_linear
    from .forcing import WindForcing
    from .grid import Grid
    from .params import Params
    from .steady import solve_steady
    cases = [
        ("Stommel", Params(delta_s=0.05), "single-gyre"),
        ("Munk free-slip", Params(delta_m=0.05), "single-gyre"),
        ("Munk no-slip", Params(delta_m=0.05, bc="no-slip", bc_ns="free-slip"), "single-gyre"),
        ("Stommel-Munk, asym. double gyre", Params(delta_s=0.02, delta_m=0.04, bc="no-slip",
                                                   bc_ns="free-slip"), "asymmetric-double-gyre"),
        ("f-plane", Params(delta_s=0.1, beta_factor=0.0), "single-gyre"),
    ]
    ns = [24, 48, 96] if args.quick else [32, 64, 128, 256]
    print("Max relative error of the finite-difference solution against the exact separable")
    print("solution, and the observed order p = log2(e_h / e_{h/2}).\n")
    print(f"{'case':<34s}" + "".join(f"{'n=' + str(n):>11s}" for n in ns) + "     orders")
    worst = math.inf
    for name, p, prof in cases:
        f = WindForcing(prof)
        errs = []
        for n in ns:
            g = Grid(n, n)
            s = solve_steady(g, p, f)
            ex = exact_linear(g.X, g.Y, p, f)["psi"]
            errs.append(float(np.max(np.abs(s.psi - ex)) / np.max(np.abs(ex))))
        orders = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
        worst = min(worst, orders[-1])
        print(f"{name:<34s}" + "".join(f"{e:>11.3e}" for e in errs) + "     "
              + " ".join(f"{o:.2f}" for o in orders))
    print(f"\nObserved order on the finest pair: >= {worst:.2f} (expected 2)")
    return 0 if worst > 1.8 else 1


def _cmd_explain(args) -> int:
    from .explain import TOPICS, explain
    if not args.topic:
        print("Topics:")
        for k, (title, _) in TOPICS.items():
            print(f"  {k:<24s}{title}")
        print("\nUsage: dewi-lanjar explain <topic>")
        return 0
    try:
        print(explain(args.topic))
    except KeyError as e:
        print(e.args[0], file=sys.stderr)
        return 2
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="dewi-lanjar",
        description="Western boundary intensification of wind-driven ocean gyres: "
                    "Stommel, Munk and inertial boundary currents with diagnostics and animations.",
        epilog="Start with: dewi-lanjar list ; dewi-lanjar run stommel_1948 ; dewi-lanjar explain overview",
    )
    ap.add_argument("--version", action="version", version=f"dewi-lanjar {__version__}")
    sub = ap.add_subparsers(dest="command")

    sp = sub.add_parser("list", help="list built-in cases")
    sp.set_defaults(func=_cmd_list)

    def add_set(p):
        p.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE",
                       help="override a case value, e.g. --set physics.delta_s=0.02 (repeatable)")

    sp = sub.add_parser("show", help="print a case and its description")
    sp.add_argument("case")
    add_set(sp)
    sp.set_defaults(func=_cmd_show)

    sp = sub.add_parser("init", help="write an editable case file from a built-in case")
    sp.add_argument("file", help="output TOML file")
    sp.add_argument("--template", default="stommel_1948", help="built-in case to start from")
    sp.add_argument("--force", action="store_true")
    add_set(sp)
    sp.set_defaults(func=_cmd_init)

    for name, fn, helptext in (("run", _cmd_run, "run one case"),
                               ("run-all", _cmd_run_all, "run all built-in cases")):
        sp = sub.add_parser(name, help=helptext)
        if name == "run":
            sp.add_argument("case", help="built-in case name or path to a TOML file")
        else:
            sp.add_argument("--kinds", default=None, help="comma list: steady,transient,sweep")
        sp.add_argument("--out", default="outputs", help="output directory (default: outputs)")
        sp.add_argument("--no-gif", action="store_true", help="skip the animation")
        sp.add_argument("--no-png", action="store_true", help="skip the static figure")
        sp.add_argument("--quiet", action="store_true")
        add_set(sp)
        sp.set_defaults(func=fn)

    sp = sub.add_parser("steady", help="solve a steady linear gyre from parameters")
    sp.add_argument("--delta-s", type=float, default=0.04, help="Stommel width / L")
    sp.add_argument("--delta-m", type=float, default=0.0, help="Munk width / L")
    sp.add_argument("--beta-factor", type=float, default=1.0, help="1 = beta plane, 0 = f plane")
    sp.add_argument("--bc", choices=["free-slip", "no-slip"], default="free-slip",
                    help="western/eastern walls")
    sp.add_argument("--bc-ns", choices=["free-slip", "no-slip"], default=None,
                    help="southern/northern walls (default: same as --bc)")
    sp.add_argument("--wind", default="single-gyre",
                    choices=["single-gyre", "double-gyre", "asymmetric-double-gyre",
                             "triple-gyre", "cyclonic-gyre"])
    sp.add_argument("--sign", type=float, default=1.0, choices=[1.0, -1.0])
    sp.add_argument("--asymmetry", type=float, default=0.5)
    sp.add_argument("--nx", type=int, default=128)
    sp.add_argument("--ny", type=int, default=None)
    sp.add_argument("--Ly", type=float, default=1.0, help="meridional extent / L")
    sp.add_argument("--png", default=None, help="also write a figure to this path")
    sp.add_argument("--why", action="store_true", help="print the derivation of the balance")
    sp.set_defaults(func=_cmd_steady)

    sp = sub.add_parser("scales", help="scales and widths from dimensional parameters")
    sp.add_argument("--L", type=float, default=5000.0, help="zonal basin width [km]")
    sp.add_argument("--lat", type=float, default=30.0, help="latitude for beta [deg]")
    sp.add_argument("--beta", type=float, default=None, help="beta [1/(m s)] (overrides --lat)")
    sp.add_argument("--r", type=float, default=1e-6, help="bottom drag [1/s]")
    sp.add_argument("--A-h", dest="A_h", type=float, default=1000.0, help="lateral viscosity [m^2/s]")
    sp.add_argument("--tau0", type=float, default=0.1, help="wind stress [N/m^2]")
    sp.add_argument("--rho0", type=float, default=1025.0, help="density [kg/m^3]")
    sp.add_argument("--H", type=float, default=1000.0, help="layer depth [m]")
    sp.set_defaults(func=_cmd_scales)

    sp = sub.add_parser("widths", help="theoretical boundary-layer widths")
    sp.add_argument("--delta-s", type=float, default=0.0)
    sp.add_argument("--delta-m", type=float, default=0.0)
    sp.add_argument("--beta-factor", type=float, default=1.0)
    sp.add_argument("--bc", choices=["free-slip", "no-slip"], default="no-slip")
    sp.set_defaults(func=_cmd_widths)

    sp = sub.add_parser("validate", help="grid-convergence test against the exact solution")
    sp.add_argument("--quick", action="store_true", help="coarser grids")
    sp.set_defaults(func=_cmd_validate)

    sp = sub.add_parser("explain", help="print a derivation")
    sp.add_argument("topic", nargs="?", default=None)
    sp.set_defaults(func=_cmd_explain)
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if not getattr(args, "command", None):
        ap.print_help()
        return 0
    try:
        return int(args.func(args) or 0)
    except (ValueError, FloatingPointError) as e:
        print(f"dewi-lanjar: error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
