# Theory

This page derives every equation that dewi-lanjar solves, starting from the
momentum equations of a rotating fluid layer. The same derivations, in plain
text, are printed by `dewi-lanjar explain <topic>`.

## 1. From momentum to vorticity

Consider a homogeneous layer of constant depth $H$ and density $\rho_0$ in a
closed basin on a $\beta$ plane, $f = f_0 + \beta y$. The flow is forced by a
wind stress $\boldsymbol{\tau} = (\tau^x, \tau^y)$ applied as a body force
over the layer, damped by linear bottom drag with coefficient $r$ and by
lateral eddy viscosity $A_h$:

$$
\begin{aligned}
\partial_t u + \mathbf{u}\cdot\nabla u - f v &= -\frac{1}{\rho_0}\partial_x p + \frac{\tau^x}{\rho_0 H} - r u + A_h \nabla^2 u,\\
\partial_t v + \mathbf{u}\cdot\nabla v + f u &= -\frac{1}{\rho_0}\partial_y p + \frac{\tau^y}{\rho_0 H} - r v + A_h \nabla^2 v,\\
\partial_x u + \partial_y v &= 0 .
\end{aligned}
$$

The continuity equation is satisfied identically by a streamfunction $\psi$,

$$
u = -\partial_y \psi, \qquad v = \partial_x \psi, \qquad
\zeta = \partial_x v - \partial_y u = \nabla^2 \psi .
$$

Taking $\partial_x$ of the $v$ equation minus $\partial_y$ of the $u$ equation
removes the pressure. The Coriolis terms give
$\partial_x(fu) + \partial_y(fv) = f\nabla\cdot\mathbf{u} + v\,\partial_y f = \beta v$,
and the advective terms give $\mathbf{u}\cdot\nabla\zeta = J(\psi, \zeta)$ with
$J(a, b) = a_x b_y - a_y b_x$. The result is the barotropic vorticity equation

$$
\partial_t \zeta + J(\psi, \zeta) + \beta\,\partial_x\psi
= \frac{\operatorname{curl}\boldsymbol{\tau}}{\rho_0 H} - r\zeta + A_h \nabla^2\zeta,
\qquad \operatorname{curl}\boldsymbol{\tau} = \partial_x\tau^y - \partial_y\tau^x .
$$

Each term has a direct meaning: advection of relative vorticity, advection of
planetary vorticity ($\beta v$), input of vorticity by the wind, removal by
bottom drag, and lateral diffusion of vorticity.

**Boundary conditions.** No flow through the walls means that each wall is a
streamline; in a simply connected basin $\psi = 0$ on all walls. Lateral
friction raises the order of the equation and needs one more condition per
wall: *free slip* (no tangential stress, $\zeta = 0$) or *no slip* (no
tangential velocity, $\partial_n\psi = 0$). dewi-lanjar lets the western and
eastern walls (`bc`) and the southern and northern walls (`bc_ns`) differ.

## 2. Scaling and the three boundary-layer widths

Scale lengths by the zonal basin width $L$, velocities by $U$ and time by
$1/(\beta L)$, and write $\operatorname{curl}\boldsymbol{\tau} = (\tau_0 / L)F$.
The terms scale as

| term | scale |
|---|---|
| $\partial_t\zeta$ | $\beta U$ |
| $J(\psi,\zeta)$ | $U^2/L^2$ |
| $\beta\,\partial_x\psi$ | $\beta U$ |
| $\operatorname{curl}\boldsymbol{\tau}/(\rho_0 H)$ | $\tau_0/(\rho_0 H L)$ |
| $r\zeta$ | $rU/L$ |
| $A_h\nabla^2\zeta$ | $A_h U/L^3$ |

The Sverdrup velocity $U = \tau_0/(\rho_0 H \beta L)$ balances the wind against
$\beta v$. Dividing by $\beta U$ gives the nondimensional equation solved by
dewi-lanjar,

$$
\boxed{\;\partial_t\zeta + \delta_I^2 J(\psi,\zeta) + b\,\partial_x\psi
= F - \delta_S\,\zeta + \delta_M^3\nabla^2\zeta,\qquad \zeta = \nabla^2\psi\;}
$$

with

$$
\delta_S = \frac{r}{\beta L},\qquad
\delta_M = \frac{1}{L}\left(\frac{A_h}{\beta}\right)^{1/3},\qquad
\delta_I = \frac{1}{L}\left(\frac{U}{\beta}\right)^{1/2},
$$

the Stommel, Munk and inertial boundary-layer widths in units of $L$. The
factor $b$ is one on a $\beta$ plane and zero in the $f$-plane control
experiment. The volume transport per unit nondimensional streamfunction is
$HUL = \tau_0/(\rho_0\beta)$, about 4.9 Sv for $\tau_0 = 0.1$ N m$^{-2}$ and
$\beta = 2\times10^{-11}$ m$^{-1}$ s$^{-1}$. The class
{class}`dewi_lanjar.DimensionalScales` and the command `dewi-lanjar scales`
perform these conversions with every number substituted.

## 3. Wind forcing

The built-in winds are zonal, $\tau^x(y)$, so $F = -\mathrm{d}\tau^x/\mathrm{d}y$.
For $\tau^x = -\cos(n\pi y/L_y)$,

$$
F(y) = -\frac{n\pi}{L_y}\sin\frac{n\pi y}{L_y} .
$$

$n = 1$ drives one anticyclonic subtropical gyre (negative curl in the
Northern Hemisphere), $n = 2$ a subtropical and a subpolar gyre, $n = 3$ three
gyres. Every built-in profile is a finite sine series
$F = \sum_n a_n \sin(k_n y)$, $k_n = n\pi/L_y$, which is what makes the exact
solution of section 7 possible. A sign of $-1$ mirrors the circulation (a
Southern Hemisphere gyre, since $\beta > 0$ in both hemispheres). Optional
periodic modulation multiplies $F$ by $1 + \epsilon\sin(2\pi t/P)$.

## 4. The Sverdrup interior

Far from the western wall the flow is broad and slow, so the time derivative,
advection and friction are small next to $b\,\partial_x\psi$. The steady
equation reduces to the Sverdrup balance

$$
b\,\partial_x\psi_I = F(x, y) .
$$

It is first order in $x$ and so admits one zonal boundary condition.
Integrating from the eastern wall, where $\psi = 0$,

$$
\psi_I(x, y) = -\frac{1}{b}\int_x^{L_x} F(x', y)\,\mathrm{d}x' .
$$

For the single-gyre wind, $\psi_I = (L_x - x)(\pi/L_y)\sin(\pi y/L_y)$:
southward flow over the whole interior with transport $\psi_I(0, y)$ at each
latitude. Because $\psi_I(0, y)\neq0$, the interior alone violates the
condition on the western wall; a thin boundary layer, in which the neglected
terms are large, must carry the interior transport back north.

## 5. Stommel's boundary layer

Keep drag only ($\delta_M = \delta_I = 0$). The steady equation is

$$
\delta_S\nabla^2\psi + b\,\partial_x\psi = F .
$$

In a layer of width $\delta_S$ the $x$ derivatives dominate:
$\partial_{xx}\psi \sim \psi/\delta_S^2 \gg \partial_{yy}\psi$, and $F$ is
small next to the two retained terms. Then
$\delta_S\,\partial_{xx}\psi + b\,\partial_x\psi = 0$, so
$\partial_x\psi \propto e^{mx}$ with $m = -b/\delta_S$. For $b > 0$ the root is
negative: the solution decays *eastward* away from a wall. An eastern-wall
layer would need $e^{m(x - L_x)}$ with $m > 0$, which the balance does not
allow. The boundary layer can only sit at the western wall. Matching to the
interior gives the composite

$$
\psi \simeq \psi_I(x, y)\left[1 - e^{-bx/\delta_S}\right],
$$

whose northward velocity is largest at the wall and decays to $1/e$ of its
peak at $x = \delta_S/b$.

**Vorticity argument.** The anticyclonic wind supplies negative vorticity
everywhere. In the interior, water moving south loses planetary vorticity at
the same rate. To return north it must gain planetary vorticity
($\beta v > 0$) and therefore needs a source of positive vorticity. Drag
contributes $-\delta_S\zeta$, positive where $\zeta < 0$. In a northward
current against the western wall $v$ decreases eastward, $\zeta = v_x < 0$,
and drag supplies the needed vorticity. Against an eastern wall $v$ would
increase eastward, $\zeta > 0$, and drag would remove vorticity of the wrong
sign. The diagnostics of section 10 show this balance directly.

## 6. Munk's boundary layer

Keep lateral friction only. The steady equation is

$$
\delta_M^3\nabla^4\psi - b\,\partial_x\psi = -F .
$$

In the layer, $\delta_M^3\,\partial_x^4\psi = b\,\partial_x\psi$, so
$\partial_x\psi\propto e^{mx}$ with $\delta_M^3 m^3 = b$:

$$
m = \frac{b^{1/3}}{\delta_M}\left\{1,\ e^{2\pi i/3},\ e^{-2\pi i/3}\right\}.
$$

The complex pair has $\operatorname{Re} m = -b^{1/3}/(2\delta_M) < 0$ and
decays eastward with oscillations of wavenumber $\sqrt3\,b^{1/3}/(2\delta_M)$.
Two decaying roots allow two wall conditions. With $\xi = x/\delta_M$ and
$b = 1$,

$$
\begin{aligned}
\text{no slip } (\psi = \psi_x = 0):&\quad
\psi = \psi_I\left[1 - e^{-\xi/2}\left(\cos\tfrac{\sqrt3\xi}{2} + \tfrac{1}{\sqrt3}\sin\tfrac{\sqrt3\xi}{2}\right)\right],\\
\text{free slip } (\psi = \psi_{xx} = 0):&\quad
\psi = \psi_I\left[1 - e^{-\xi/2}\left(\cos\tfrac{\sqrt3\xi}{2} - \tfrac{1}{\sqrt3}\sin\tfrac{\sqrt3\xi}{2}\right)\right].
\end{aligned}
$$

Differentiating, the no-slip velocity is
$v\propto e^{-\xi/2}\sin(\sqrt3\xi/2)$, which first changes sign at
$x = 2\pi\delta_M/\sqrt3 \approx 3.63\,\delta_M$; the free-slip velocity
$v\propto e^{-\xi/2}[\cos(\sqrt3\xi/2) + \sin(\sqrt3\xi/2)/\sqrt3]$ first
changes sign at $x = 4\pi\delta_M/(3\sqrt3)\approx 2.42\,\delta_M$. Offshore
of that point flows a weak counter-current, a feature the Stommel layer does
not have.

**General local layer.** With drag and viscosity together the layer equation
for $X = \partial_x\psi$ is

$$
\delta_M^3 X''' - \delta_S X' - b X = 0,\qquad
\delta_M^3 m^3 - \delta_S m - b = 0 .
$$

{func}`dewi_lanjar.boundary_layer_profile` solves this problem with the wall
conditions above, and {func}`dewi_lanjar.theory_widths` measures its e-folding
width and counter-current position with the same algorithm that is applied to
the numerical solutions, so theory and model are compared like for like.

## 7. Exact separable solution of the steady linear problem

For $F = \sum_n a_n\sin(k_n y)$ and walls at $y = 0, L_y$ that are free slip
(or any walls when $\delta_M = 0$), $\sin(k_n y)$ satisfies all conditions on
those walls. Writing $\psi = \sum_n X_n(x)\sin(k_n y)$ turns the steady linear
equation into one ordinary differential equation per mode:

$$
\delta_M^3\left(D^2 - k^2\right)^2 X - \delta_S\left(D^2 - k^2\right)X - b\,DX = -a,
\qquad D = \frac{\mathrm{d}}{\mathrm{d}x}.
$$

A constant particular solution is
$X_p = -a/(\delta_M^3 k^4 + \delta_S k^2)$. Homogeneous solutions $e^{mx}$
satisfy the characteristic quartic

$$
\delta_M^3 m^4 - \left(2\delta_M^3 k^2 + \delta_S\right) m^2 - b\,m
+ \left(\delta_M^3 k^4 + \delta_S k^2\right) = 0,
$$

which reduces to the quadratic $\delta_S m^2 + b\,m - \delta_S k^2 = 0$ when
$\delta_M = 0$ (Stommel's 1948 solution). The constants follow from
$X(0) = X(L_x) = 0$ plus, when $\delta_M > 0$, $X'' = 0$ (free slip) or
$X' = 0$ (no slip) at both walls. To avoid overflow each exponential is
written as $e^{m(x - x_m)}$ with $x_m = L_x$ when $\operatorname{Re} m > 0$
and $x_m = 0$ otherwise, so that no basis function exceeds one inside the
basin; the small linear system is solved after column scaling.
{func}`dewi_lanjar.exact_linear` evaluates $\psi$, $u$, $v$ and $\zeta$ to
machine precision. The test suite checks it against an independent
50-digit evaluation of Stommel's closed form.

## 8. Inertial boundary currents

When $\delta_I$ is comparable with or larger than the frictional widths, the
advection of relative vorticity enters the layer balance. In the inviscid,
steady limit the equation states that potential vorticity
$q = \delta_I^2\zeta + b y$ is constant along streamlines,

$$
J\!\left(\psi,\ \delta_I^2\zeta + b y\right) = 0 .
$$

A fluid column carried north by the boundary current gains planetary
vorticity and must lose relative vorticity to conserve $q$. The current cannot
turn east where the linear current would: it overshoots north and closes in a
recirculation cell, so the gyre center moves north and west and the
boundary-current transport exceeds the Sverdrup value. Useful ratios,
reported by {meth}`dewi_lanjar.Params.regime`, are the boundary-layer
Reynolds numbers $(\delta_I/\delta_M)^3$ and $(\delta_I/\delta_S)^2$.

With free-slip walls and weak friction the inertial solution can grow into a
basin-filling mode close to Fofonoff's free solution, in which $q$ is nearly a
linear function of $\psi$; the kinetic energy then rises far above its linear
value and the western intensification is lost (the case `inertial_runaway`).
No-slip walls generate vorticity of the opposite sign at the wall and arrest
this runaway (the case `inertial_gyre`). In a double-gyre basin the two
boundary currents meet and separate as a free jet that becomes unstable and
sheds eddies (the case `double_gyre_eddies`).

## 9. Spin-up by Rossby waves

Linearize about rest and drop friction:
$\partial_t\nabla^2\psi + b\,\partial_x\psi = F$. Plane waves
$e^{i(kx + ly - \omega t)}$ obey the Rossby dispersion relation

$$
\omega = -\frac{b\,k}{k^2 + l^2} .
$$

Long waves ($k^2 \ll l^2$) are nondispersive with zonal group velocity
$-b/l^2$: westward. After the wind is switched on, the Sverdrup interior is
established behind a front of long Rossby waves that leaves the eastern wall
and crosses the basin in a time of order $L_x l^2/b$ (in units of
$1/(\beta L)$). Long waves reflected at the western wall become short waves
with small eastward group velocity, so their energy accumulates near the
western wall, where friction finally balances it. The case
`stommel_munk_spinup` shows this process; `pulsating_wind` drives it
periodically.

## 10. Energy and vorticity budgets

**Energy.** Multiply the vorticity equation by $-\psi$ and integrate over the
basin. With $\psi = 0$ on the walls, $-\int\psi\,\partial_t\zeta = \frac{\mathrm{d}}{\mathrm{d}t}\frac12\int|\nabla\psi|^2$,
$\int\psi J(\psi,\zeta) = 0$ and $\int\psi\,\partial_x\psi = 0$, so

$$
\frac{\mathrm{d}KE}{\mathrm{d}t} = W - 2\delta_S KE - \delta_M^3\int\zeta^2\,\mathrm{d}A,
\qquad W = -\int\psi F\,\mathrm{d}A .
$$

The boundary terms of the viscous integral vanish for both free-slip
($\zeta = 0$) and no-slip ($\partial_n\psi = 0$) walls. In a steady state the
wind work equals the dissipation.

**Vorticity.** Integrating the equation over the basin, the beta term gives
$b\int[\psi(L_x, y) - \psi(0, y)]\,\mathrm{d}y = 0$: the planetary vorticity
gradient redistributes vorticity but cannot create or destroy it. In a steady
state the net vorticity input of the wind must therefore be removed by
friction. Splitting a gyre's latitude band into a western strip and the rest
separates the two balances:

$$
\text{interior:}\ \ F \approx b\,\partial_x\psi,\qquad
\text{western strip:}\ \ b\,\partial_x\psi \approx -\delta_S\zeta + \delta_M^3\nabla^2\zeta .
$$

On the grid the beta term summed along a row equals
$b(\psi_{N-1} - \psi_1)/2$ rather than zero, a truncation error of order $h$
at the walls; dewi-lanjar reports it rather than hiding it.

## 11. The f-plane control

With $b = 0$ the steady drag-only problem is $\delta_S\nabla^2\psi = F$. For a
wind that does not depend on $x$ the solution is symmetric about the middle of
the basin: no Sverdrup interior, no boundary current, an intensification index
of one and a kinetic-energy centroid at $x = L_x/2$. This experiment isolates
$\beta$ as the cause of the asymmetry; the sweep case `beta_emergence` shows
the transition continuously as $b$ grows from zero to one.

## 12. Concentration of kinetic energy

Let $e_{ij} = \frac12(u^2 + v^2)$ at the grid nodes and $p_{ij} = e_{ij}/\sum e$.
The normalized Shannon entropy

$$
H = -\frac{1}{\ln N}\sum_{ij} p_{ij}\ln p_{ij}
$$

equals one for kinetic energy spread uniformly and decreases as energy
concentrates. The same quantity for the zonal marginal $\sum_j e_{ij}$
isolates east-west concentration, and the zonal centroid
$x_c = \sum_i x_i\sum_j e_{ij} / (L_x\sum e)$ moves from 0.5 toward zero as a
western boundary current forms. These numbers summarize in one scalar what a
map shows: western intensification is a concentration of kinetic energy, and
its strength can be tracked through parameter sweeps and spin-up.
