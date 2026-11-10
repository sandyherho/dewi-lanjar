# Diagnostics

{func}`dewi_lanjar.diagnostics.compute` returns a dictionary of physically
interpretable numbers for any state; {func}`dewi_lanjar.diagnostics.describe`
formats it. Every case run writes the same numbers to `summary.txt` and
`diagnostics.json`, and transient runs record the scalar ones as time series.

All integrals are sums over interior nodes times $\Delta x\,\Delta y$, the
inner product under which the discrete Laplacian is symmetric. With this
choice the energy and vorticity budgets close to round-off, so a residual
measures departure from a steady state, not quadrature error.

## Energetics

`KE`, `enstrophy`
: $KE = -\frac12\sum\psi\zeta\,\Delta A$, which equals the discrete
  $\frac12\int|\nabla\psi|^2$, and $Z = \frac12\sum\zeta^2\,\Delta A$.

`energy_wind`, `energy_drag`, `energy_viscous`, `energy_advection`, `energy_beta`, `dKE_dt`, `energy_imbalance`
: The terms of $\mathrm dKE/\mathrm dt = W + D_{\text{drag}} + D_{\text{visc}} + A_{\text{adv}} + A_\beta$,
  each computed as $-\sum\psi\,(\text{term})\,\Delta A$. The advective and beta
  terms vanish identically for the Arakawa Jacobian and the centered
  difference and are reported to show that they do. `energy_imbalance` is
  $\mathrm dKE/\mathrm dt$ divided by the wind work: zero to round-off in a
  steady state.

## Vorticity budget

`vorticity_budget`
: The integrals of wind, beta, drag, viscosity and advection over four
  regions: the whole basin, the latitude band of the diagnosed gyre
  (`gyre_band`), the western strip of that band (`strip_width`, by default
  three boundary-layer decay scales) and the rest of the band. In a steady
  linear state the interior shows the Sverdrup balance (wind against beta)
  and the western strip the boundary-layer balance (beta against drag and
  viscosity). For a double gyre the band stops at the inter-gyre boundary, so
  the two gyres do not cancel.

## Boundary current

The section is taken at `section_y`, the latitude of the largest $|\psi|$.

`wbc_v_peak`, `wbc_x_peak`
: Peak along-wall velocity and its distance from the western wall (zero for
  free slip, inside the domain for no slip).

`wbc_e_fold`, `theory_e_fold`
: Distance from the wall at which $|v|$ has fallen to $|v_{\text{peak}}|/e$
  east of the peak, for the model and for the boundary-layer theory of
  {func}`dewi_lanjar.boundary_layer_profile`. For a pure Stommel layer the
  theoretical value is exactly $\delta_S$.

`wbc_zero_crossing`, `theory_zero_crossing`
: First sign change of $v$: the inner edge of the offshore counter-current.
  Munk theory gives $2\pi\delta_M/\sqrt3$ (no slip) and
  $4\pi\delta_M/(3\sqrt3)$ (free slip); a Stommel layer has none.

`wbc_transport`, `sverdrup_transport`, `transport_ratio`
: The streamfunction at the outer edge of the current (its volume
  transport), the Sverdrup prediction $\psi_I$ at the same point, and their
  ratio. Below one: interior friction is not negligible (large $\delta$).
  Above one: a recirculation adds transport (inertial regime). The ratio is
  reported as `nan` where $\psi_I$ nearly vanishes.

## Asymmetry and concentration

`intensification_index`
: Peak $|v|$ in the western half of the section over peak $|v|$ in the
  eastern half. One for a symmetric gyre; above ten for typical western
  boundary currents.

`ke_centroid_x`
: Zonal center of mass of kinetic energy in units of $L_x$; 0.5 when
  symmetric.

`ke_entropy`, `ke_zonal_entropy`
: Normalized Shannon entropies of the kinetic-energy distribution over the
  basin and of its zonal marginal ([theory, section 12](theory.md)). One for
  uniform energy; they decrease as energy concentrates in the boundary layer.

## Resolution

`points_per_layer`
: Grid points per boundary-layer decay scale $1/|\operatorname{Re} m|$.
  Below three the summary warns.

`grid_scale_fraction`
: Fraction of the vorticity norm removed by a $3\times3$ running mean. Above
  0.3 the summary warns that the smallest scales are dominated by
  two-grid-interval noise.

## Regime

`regime`
: The output of {meth}`dewi_lanjar.Params.regime`: the three widths, the
  ratios $\delta_I/\delta_S$ and $\delta_I/\delta_M$, the boundary-layer
  Reynolds numbers and a label (linear Stommel, linear Munk, Stommel-Munk,
  weakly nonlinear, inertial, f-plane).
