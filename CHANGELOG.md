# Changelog

## 0.0.1 (2026-11-09)

First release.

* Barotropic vorticity equation on a beta plane in a closed rectangular basin,
  nondimensionalized by the Stommel, Munk and inertial widths.
* Steady linear sparse direct solver; time-dependent nonlinear solver
  (Arakawa Jacobian, sine-transform Poisson inversion, SSP-RK3, adaptive time
  step); free-slip and no-slip walls, set separately for the meridional and
  zonal walls.
* Exact separable solution of the steady linear problem for any sine-series
  wind; Sverdrup interior; boundary-layer theory for any combination of drag
  and viscosity.
* Diagnostics: energy and vorticity budgets, boundary-current widths against
  theory, transport closure, intensification index, kinetic-energy centroid
  and Shannon entropies, resolution and grid-noise checks.
* Lagrangian floats.
* Twelve built-in cases (steady, transient, sweep) in TOML, GIF animations and
  figures, command-line interface with `run`, `steady`, `scales`, `widths`,
  `validate` and `explain`.
* Sphinx documentation for Read the Docs with full derivations.
