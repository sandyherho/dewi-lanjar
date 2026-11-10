# Quick start

## From the command line

```bash
dewi-lanjar list                       # the built-in cases
dewi-lanjar run stommel_1948           # writes outputs/stommel_1948/
dewi-lanjar explain stommel            # the derivation behind it
```

Each run writes to `outputs/<case>/`:

| file | content |
|---|---|
| `<case>.gif` | the animation |
| `<case>.png` | static figure: streamfunction, vorticity, velocity section, budget |
| `summary.txt` | case description and every diagnostic, in physical units when available |
| `diagnostics.json` | the same numbers, machine readable |
| `fields.npz` | final fields, coordinates and (transient runs) time series |
| `case.toml` | the fully resolved configuration, to reproduce or edit the run |

Change any parameter without writing a file:

```bash
dewi-lanjar run stommel_1948 --set physics.from_scales=false --set physics.delta_s=0.02
dewi-lanjar run munk_1950 --set physics.bc=free-slip --set grid.nx=256 --set grid.ny=256
dewi-lanjar run inertial_gyre --set physics.delta_i=0.05 --out runs
```

Start your own case from a built-in one:

```bash
dewi-lanjar init my_gyre.toml --template double_gyre_eddies
# edit my_gyre.toml
dewi-lanjar run my_gyre.toml
```

Quick steady solutions and ocean scales:

```bash
dewi-lanjar steady --delta-s 0.02 --delta-m 0.03 --bc no-slip --bc-ns free-slip --why
dewi-lanjar scales --L 6000 --lat 35 --r 1e-6 --A-h 2000
dewi-lanjar widths --delta-m 0.03 --bc no-slip
```

## From Python

```python
from dewi_lanjar import Grid, Params, WindForcing, solve_steady, exact_linear, diagnostics

grid = Grid(192, 192)
p = Params(delta_s=0.02, delta_m=0.03, bc="no-slip", bc_ns="free-slip")
wind = WindForcing("double-gyre")

sol = solve_steady(grid, p, wind)                         # sparse direct solve
exact = exact_linear(grid.X, grid.Y, p, wind)["psi"]      # exact separable solution
print(abs(sol.psi - exact).max() / abs(exact).max())      # about 1e-3, second order

d = diagnostics.compute(sol.psi, sol.zeta, grid, p, wind)
print(diagnostics.describe(d))
```

Time-dependent, nonlinear runs:

```python
from dewi_lanjar import TransientSolver, RunConfig

p = Params(delta_s=0.02, delta_m=0.02, delta_i=0.035, bc="no-slip")
run = TransientSolver(Grid(192, 192), p, WindForcing()).run(
    RunConfig(t_end=200, n_frames=100, n_floats=500)
)
run.series["KE"]            # kinetic-energy time series
run.final["transport_ratio"]
run.psi.shape               # (100, 193, 193) stored snapshots
```

Dimensional parameters:

```python
from dewi_lanjar import DimensionalScales

ocean = DimensionalScales(L=5e6, beta=2e-11, r=1e-6, A_h=1000, tau0=0.1, H=1000)
print(ocean.explain())       # every scale with numbers substituted
p = ocean.to_params()
ocean.sverdrups(3.14)        # nondimensional psi -> Sv
```

Cases from Python:

```python
from dewi_lanjar import load_case, run_case

case = load_case("beta_emergence", ["sweep.n=30"])
result = run_case(case, "outputs")
result.files["gif"]
```
