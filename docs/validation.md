# Validation

Every number in dewi-lanjar can be checked. This page lists what is verified,
how, and the results on the reference machine. All of it is reproduced by
`dewi-lanjar validate`, `python tools/validation_report.py` and `pytest`.

## Grid convergence against the exact solution

The finite-difference steady solution is compared with the exact separable
solution ([theory, section 7](theory.md)) on grids of $n\times n$ cells. The
table gives the maximum error relative to $\max|\psi|$ and the observed order
$p = \log_2(e_h/e_{h/2})$.

| case | $n=32$ | $n=64$ | $n=128$ | $n=256$ | orders |
|---|---|---|---|---|---|
| Stommel, $\delta_S = 0.05$ | 1.42e-2 | 3.57e-3 | 8.85e-4 | 2.21e-4 | 2.00 2.01 2.00 |
| Munk free slip, $\delta_M = 0.05$ | 1.49e-2 | 3.61e-3 | 8.95e-4 | 2.23e-4 | 2.05 2.01 2.00 |
| Munk no slip, $\delta_M = 0.05$ | 3.73e-2 | 9.39e-3 | 2.35e-3 | 5.89e-4 | 1.99 2.00 2.00 |
| Stommel-Munk, asymmetric double gyre | 4.72e-2 | 1.32e-2 | 3.35e-3 | 8.40e-4 | 1.84 1.98 1.99 |
| f plane, $\delta_S = 0.1$ | 5.21e-5 | 1.31e-5 | 3.28e-6 | 8.20e-7 | 1.99 1.99 2.00 |

The order is two in every configuration, including no-slip walls with the
Thom closure.

## What the test suite checks

* **Exact solution.** The separable solution agrees with an independent
  50-digit evaluation of Stommel's closed form to a relative $10^{-11}$;
  property-based tests (Hypothesis) check that random modes satisfy the
  ordinary differential equation and every wall condition, for random widths,
  $\beta$ factors and boundary conditions.
* **Boundary-layer theory.** The numerical theory widths reproduce the closed
  forms $\delta_S$, $2\pi\delta_M/\sqrt3$ and $4\pi\delta_M/(3\sqrt3)$, and the
  numerical profile reproduces Munk's formula to $10^{-12}$.
* **Operators.** The Arakawa Jacobian conserves discrete energy (and
  enstrophy, when $\zeta$ vanishes on the walls) to round-off on random
  fields and grids, is antisymmetric, and is second-order accurate; the sine
  transform Poisson solver inverts the five-point Laplacian to round-off.
* **Solver consistency.** A long linear time integration converges to the
  sparse steady solution to $10^{-8}$, for free-slip and no-slip walls.
* **Conservation.** Without forcing and friction, kinetic energy is
  conserved to $10^{-4}$ relative over the test interval (the residual is the
  time-stepping error); during a forced spin-up the time integral of the
  energy budget matches the change of kinetic energy.
* **Physics.** The f-plane gyre is exactly symmetric (centroid 0.5,
  intensification index 1) while the beta-plane gyre is not; mirroring the
  wind mirrors the solution; in a steady state the energy budget closes to
  $10^{-9}$ and the interior and western-strip vorticity balances have the
  expected structure; measured Munk widths match theory within a few
  percent; inertia moves the gyre north and raises its transport.
* **Floats.** In a steady flow floats stay on their streamlines.
* **Interface.** Every built-in case loads and validates, case files round
  trip through TOML, every command of the CLI runs, small steady, transient
  and sweep cases write all their outputs.

## Limitations

* Second-order accuracy needs resolution: keep three or more grid points per
  boundary-layer decay scale and watch the grid-scale vorticity fraction.
* The grid is uniform; very thin layers ($\delta \lesssim 0.005$) need many
  points.
* The model is a single homogeneous layer. Baroclinic effects, topography,
  coastline geometry and the real wind field are outside its scope; the
  numbers are those of the idealized problem, not of a particular ocean.
* The exact separable solution requires free-slip southern and northern
  walls when lateral friction is present.
* The steady solver is linear. Nonlinear states are obtained by time
  integration.
