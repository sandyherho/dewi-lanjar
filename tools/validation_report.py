"""Reproduce the numbers on the Validation page of the documentation.

Usage: python tools/validation_report.py
"""

import math

import numpy as np

from dewi_lanjar import (
    Grid, Params, RunConfig, TransientSolver, WindForcing, diagnostics, exact_linear,
    solve_steady, theory_widths,
)
from dewi_lanjar.cli import main as cli


def section(title):
    print("\n" + title + "\n" + "-" * len(title))


section("1. Grid convergence against the exact separable solution")
cli(["validate"])

section("2. Boundary-layer theory against closed forms")
for dm in (0.02, 0.05):
    w = theory_widths(Params(delta_m=dm, bc="no-slip"))
    print(f"Munk no slip  dM={dm}: counter-current {w['zero_crossing']:.10f}  "
          f"closed form {2 * math.pi * dm / math.sqrt(3):.10f}")
    w = theory_widths(Params(delta_m=dm, bc="free-slip"))
    print(f"Munk free slip dM={dm}: counter-current {w['zero_crossing']:.10f}  "
          f"closed form {4 * math.pi * dm / (3 * math.sqrt(3)):.10f}")
w = theory_widths(Params(delta_s=0.04))
print(f"Stommel dS=0.04: e-fold {w['e_fold']:.10f}  closed form 0.04")

section("3. Measured widths of numerical solutions (256 x 64 grid)")
for p in (Params(delta_s=0.03), Params(delta_m=0.03, bc="no-slip", bc_ns="free-slip"),
          Params(delta_m=0.03, bc="free-slip")):
    g = Grid(256, 64)
    s = solve_steady(g, p, WindForcing())
    d = diagnostics.compute(s.psi, s.zeta, g, p, s.forcing)
    print(f"{p.regime()['label']:<28s} bc={p.bc:<9s} e-fold {d['wbc_e_fold']:.4f} "
          f"(theory {d['theory_e_fold']:.4f})  counter-current {d['wbc_zero_crossing']:.4f} "
          f"(theory {d['theory_zero_crossing']:.4f})")

section("4. Steady-state budgets (Stommel-Munk, double gyre, 128 x 128)")
g = Grid(128, 128)
p = Params(delta_s=0.02, delta_m=0.02, bc="no-slip")
f = WindForcing("double-gyre")
s = solve_steady(g, p, f)
d = diagnostics.compute(s.psi, s.zeta, g, p, f)
print(f"energy imbalance (dKE/dt / W): {d['energy_imbalance']:.2e}")
for r in ("west", "interior"):
    vb = d["vorticity_budget"][r]
    print(f"{r:<9s} wind {vb['wind']: .4f}  beta {vb['beta']: .4f}  drag {vb['drag']: .4f}  "
          f"viscosity {vb['viscosity']: .4f}  sum {vb['tendency']: .1e}")

section("5. Linear time integration converges to the sparse solution")
g = Grid(48, 48)
p = Params(delta_s=0.08, delta_m=0.04, bc="no-slip")
st = solve_steady(g, p, f)
r = TransientSolver(g, p, f).run(RunConfig(t_end=600, n_frames=2, store_dtype="float64"))
print(f"max |psi_transient - psi_steady| / max|psi| = "
      f"{np.max(np.abs(r.psi[-1] - st.psi)) / np.max(np.abs(st.psi)):.2e}")

section("6. f plane versus beta plane (Stommel, 128 x 128)")
for b in (0.0, 1.0):
    p = Params(delta_s=0.04, beta_factor=b)
    s = solve_steady(Grid(128, 128), p, WindForcing())
    d = diagnostics.compute(s.psi, s.zeta, s.grid, p, s.forcing)
    print(f"b={b}: intensification {d['intensification_index']:.3f}  centroid "
          f"{d['ke_centroid_x']:.3f}  zonal KE entropy {d['ke_zonal_entropy']:.3f}")

ex = exact_linear(np.array(0.5), np.array(0.5), Params(delta_s=0.04), WindForcing())["psi"]
print(f"\nexact psi at basin center (Stommel, dS=0.04): {float(ex):.12f}")
