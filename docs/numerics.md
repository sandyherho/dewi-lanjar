# Numerical methods

dewi-lanjar uses second-order finite differences on a uniform grid,
an exact fast Poisson solver, an energy-conserving Jacobian and a
strong-stability-preserving Runge-Kutta scheme. The choices were made so
that each discrete operator inherits the property of its continuous
counterpart that matters for this problem, and so that every result can be
checked against an exact reference.

## Grid

The basin $[0, L_x]\times[0, L_y]$ carries $(n_x + 1)\times(n_y + 1)$ nodes
$x_i = i\,\Delta x$, $y_j = j\,\Delta y$. Arrays are indexed `[i, j]` with `i`
running west to east. The streamfunction is zero on the walls, so the
$(n_x - 1)\times(n_y - 1)$ interior nodes carry the unknowns.

## Spatial operators

Five-point Laplacian and centered first differences:

$$
(\nabla_h^2 f)_{ij} = \frac{f_{i+1,j} - 2f_{ij} + f_{i-1,j}}{\Delta x^2}
+ \frac{f_{i,j+1} - 2f_{ij} + f_{i,j-1}}{\Delta y^2},\qquad
(\delta_x f)_{ij} = \frac{f_{i+1,j} - f_{i-1,j}}{2\Delta x}.
$$

**Arakawa Jacobian.** Advection uses the nine-point Jacobian of Arakawa
(1966), the average $J_A = \frac13(J^{++} + J^{+\times} + J^{\times+})$ of three
second-order forms (written out in {func}`dewi_lanjar.operators.arakawa_jacobian`).
With $\psi$ constant on the boundary it satisfies
$\sum_{ij}\psi_{ij}J_A(\psi,\zeta)_{ij} = 0$ exactly, so advection neither
creates nor destroys the discrete kinetic energy $-\frac12\sum\psi\zeta$; for
$\zeta$ also vanishing on the boundary it conserves enstrophy too. This
removes the nonlinear (aliasing) instability that makes naive Jacobians blow
up in long integrations. The test suite verifies both identities to round-off
on random fields.

**Wall vorticity.** Along the western wall $\psi = 0$ implies
$\partial_{yy}\psi = 0$, so $\zeta = \partial_{xx}\psi$ there.

* Free slip: $\zeta_{0j} = 0$.
* No slip: $\partial_x\psi = 0$ is imposed with the ghost value
  $\psi_{-1,j} = \psi_{1,j}$, which gives Thom's (1933) formula
  $\zeta_{0j} = 2\psi_{1j}/\Delta x^2$.

The other walls are analogous. The conditions on the western and eastern
walls (`bc`) and on the southern and northern walls (`bc_ns`) can differ.

## Poisson inversion

$\nabla_h^2\psi = \zeta$ with $\psi = 0$ on the walls is solved exactly (to
round-off) by the type-I discrete sine transform. The vectors
$\sin(p\pi i/n_x)\sin(q\pi j/n_y)$ are eigenvectors of the Dirichlet
five-point Laplacian with eigenvalues

$$
\lambda_{pq} = \frac{2\cos(p\pi/n_x) - 2}{\Delta x^2} + \frac{2\cos(q\pi/n_y) - 2}{\Delta y^2},
$$

so $\psi = \mathrm{DST}^{-1}\left[\mathrm{DST}[\zeta]/\lambda\right]$ in
$O(N\log N)$ operations.

## Steady linear solver

The steady linear equation becomes the sparse system

$$
\left[\delta_M^3\left(L L + B\right) - \delta_S L - b\,D_x\right]\boldsymbol\psi = -\mathbf F ,
$$

where $L$ is the Dirichlet Laplacian, $D_x$ the centered difference and $B$ a
diagonal matrix carrying the Thom wall vorticity into the stencil of the first
interior row and column ($B = 0$ on free-slip walls). It is assembled with
Kronecker products and solved with SuperLU. These are the operators of the
time stepper, so a linear time integration relaxes to this solution to
round-off (verified in the tests to $10^{-8}$).

## Time integration

The three-stage, third-order strong-stability-preserving Runge-Kutta method
(Shu and Osher 1988),

$$
\begin{aligned}
\zeta^{(1)} &= \zeta^n + \Delta t\,R(\zeta^n),\\
\zeta^{(2)} &= \tfrac34\zeta^n + \tfrac14\left[\zeta^{(1)} + \Delta t\,R(\zeta^{(1)})\right],\\
\zeta^{n+1} &= \tfrac13\zeta^n + \tfrac23\left[\zeta^{(2)} + \Delta t\,R(\zeta^{(2)})\right],
\end{aligned}
$$

where every evaluation of $R$ inverts the Poisson problem and refreshes the
wall vorticity. Its stability region contains the imaginary axis up to
$\sqrt3$ and the negative real axis up to about 2.51. The time step is the
smallest of

| process | bound |
|---|---|
| advection | $\sqrt3 / \left[\delta_I^2\max(\lvert u\rvert/\Delta x + \lvert v\rvert/\Delta y)\right]$ |
| Rossby waves | $\sqrt3 / (b L_y / 2\pi)$, since $\lvert\omega\rvert = b\lvert k_x\rvert/k^2 \le b/(2 k_y^{\min})$ |
| viscosity | $\frac12\cdot 2.51 / \left[\delta_M^3(4/\Delta x^2 + 4/\Delta y^2)\right]$ |
| drag | $2.51/\delta_S$ |

multiplied by the safety factor `cfl`. The step is shortened when needed so
that output times are hit exactly.

## Lagrangian floats

Floats move with $\mathrm{d}\mathbf x/\mathrm dt = \delta_I^2\mathbf u$
(the velocity unit times the time unit, in units of $L$). Node velocities come
from second-order differences (one-sided at the walls), are interpolated
bilinearly and integrated with Heun's method; in time-dependent runs the
predictor uses the old and the corrector the new velocity field.

## Initial perturbation

A symmetric wind and a symmetric initial state remain exactly symmetric under
a symmetric discretization. `run.init_noise` adds a smooth random initial
vorticity field (seeded, so runs are reproducible) that breaks such
symmetries, which is needed for the instability of a double-gyre jet to
develop.

## Resolution

The second-order error is controlled by the number of grid points across the
narrowest feature. The diagnostics report

* the number of grid points per boundary-layer decay scale, with a warning
  below three; and
* the grid-scale fraction of the vorticity (the part removed by a
  $3\times3$ running mean), with a warning above 0.3.

Stability alone is not accuracy: thanks to the Arakawa Jacobian an
under-resolved eddying run stays bounded, but it fills with
two-grid-interval noise. The second diagnostic exposes that case.
