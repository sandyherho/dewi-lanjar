"""Plain-text derivations printed by ``dewi-lanjar explain <topic>``.

The documentation contains the same derivations typeset in LaTeX; these
versions use ASCII notation so they read well in a terminal.
"""

from __future__ import annotations

TOPICS: dict[str, tuple[str, str]] = {}


def _topic(name: str, title: str):
    def deco(fn):
        TOPICS[name] = (title, fn())
        return fn
    return deco


@_topic("overview", "What dewi-lanjar solves")
def _overview():
    return """
Wind blowing over a closed ocean basin drives a gyre. The gyre is not
symmetric: on the western side of every ocean basin the return flow is
squeezed into a narrow, fast boundary current (Gulf Stream, Kuroshio,
Agulhas, Brazil, East Australian Current), while the eastern side has a
broad, slow drift. dewi-lanjar explains and computes this asymmetry.

Model: a single homogeneous layer of depth H on a beta plane,
f = f0 + beta*y, in a rectangular basin of zonal width L. Taking the curl of
the momentum equations eliminates pressure and leaves the barotropic
vorticity equation for the streamfunction psi (u = -psi_y, v = psi_x):

    d(zeta)/dt + J(psi, zeta) + beta * v
        = curl(tau)/(rho0 H) - r * zeta + A_h * lap(zeta),   zeta = lap(psi)

    J(a, b) = a_x b_y - a_y b_x           (advection of relative vorticity)
    beta * v                              (advection of planetary vorticity)
    curl(tau)/(rho0 H)                    (wind input of vorticity)
    -r * zeta                             (linear bottom drag)
    A_h * lap(zeta)                       (lateral eddy viscosity)

Topics: nondimensionalization, sverdrup, stommel, munk, inertial,
fplane, budgets, numerics, diagnostics, floats.
"""


@_topic("nondimensionalization", "Scaling and the three boundary-layer widths")
def _nondim():
    return """
Scale lengths by L, velocities by U, time by T = 1/(beta L):

    x = L x',   psi = U L psi',   zeta = (U/L) zeta',   t = t'/(beta L),
    curl(tau) = (tau0/L) F'.

Insert into the vorticity equation term by term:

    d(zeta)/dt        ~ (U/L)(beta L)       = beta U
    J(psi, zeta)      ~ U^2 / L^2
    beta v            ~ beta U
    curl(tau)/(rho0H) ~ tau0/(rho0 H L)
    r zeta            ~ r U / L
    A_h lap(zeta)     ~ A_h U / L^3

Choose U so that the wind is balanced by beta*v (the Sverdrup scale):

    tau0/(rho0 H L) = beta U   =>   U = tau0 / (rho0 H beta L).

Divide every term by beta U:

    zeta_t + dI^2 J(psi, zeta) + b psi_x = F - dS zeta + dM^3 lap(zeta)

    dS = r / (beta L)               Stommel width / L
    dM = (A_h / beta)^(1/3) / L     Munk width / L
    dI = (U / beta)^(1/2) / L       inertial width / L
    b  = 1 (beta plane) or 0 (f-plane control experiment)

Each width is the distance over which beta*v is balanced by the term that
carries it: the boundary layer is thin when friction and inertia are weak.
Run `dewi-lanjar scales` to evaluate these numbers for an ocean.
"""


@_topic("sverdrup", "Sverdrup balance in the interior")
def _sverdrup():
    return """
Far from the western wall the flow is slow and broad, so friction,
advection and the time derivative are small compared with beta*v. The
steady equation reduces to the Sverdrup (1947) balance

    b psi_x = F(x, y)      (dimensional: beta v = curl(tau)/(rho0 H)).

It is a first-order ODE in x, so it admits ONE zonal boundary condition.
Integrating from the eastern wall, where psi = 0,

    psi_I(x, y) = -(1/b) * integral_x^Lx F(x', y) dx'.

For the single-gyre wind F = -(pi/Ly) sin(pi y/Ly) this gives
psi_I = (Lx - x)(pi/Ly) sin(pi y/Ly) > 0: southward flow everywhere in
the interior. psi_I does not vanish at x = 0, so the mass carried south by
the interior has nowhere to go. A thin western boundary layer, where the
neglected terms matter, must carry it back north. The Sverdrup transport
psi_I(0, y) is the transport that boundary current has to carry.

Why the WEST? See `dewi-lanjar explain stommel` for the boundary-layer root
and the vorticity argument.
"""


@_topic("stommel", "Stommel (1948): bottom drag")
def _stommel():
    return """
Steady, linear, drag only (dM = dI = 0):

    dS lap(psi) + b psi_x = F.

Boundary-layer argument. Near x = 0 the layer coordinate xi = x/dS
stretches the x derivatives: psi_xx ~ psi/dS^2 >> psi_yy and F is
negligible compared with the two remaining terms. Then

    dS psi_xx + b psi_x = 0   =>   psi_x ~ exp(m x),  dS m + b = 0,  m = -b/dS.

With b > 0, m < 0: the solution DECAYS eastward away from the wall. A layer
at the eastern wall would need psi_x ~ exp(m (x - Lx)) decaying westward,
i.e. m > 0, which the balance does not allow. The boundary layer can only
sit on the west. Matching to the interior:

    psi(x, y) ~ psi_I(x, y) * [1 - exp(-b x / dS)].

The northward velocity v = psi_x is largest at the wall and falls to 1/e of
its peak at x = dS/b: the e-folding width is exactly dS.

Vorticity view. The anticyclonic wind supplies negative vorticity
everywhere. In the interior, water moving south loses planetary vorticity
at the same rate (Sverdrup). To return north, water must GAIN planetary
vorticity (beta*v > 0), so it needs a source of positive vorticity. Drag
supplies -dS*zeta, which is positive only where zeta < 0. In a northward
current against the WESTERN wall v decreases eastward, so zeta = v_x < 0
and drag supplies exactly the needed positive vorticity. Against an eastern
wall v would increase eastward, zeta > 0, and drag would remove vorticity
of the wrong sign. Integrated over the gyre band, the wind input is taken
out by drag in the western layer (see `dewi-lanjar explain budgets`).

Exact solution. For F = sum a_n sin(k_n y), psi = sum X_n(x) sin(k_n y) with

    dS (X'' - k^2 X) + b X' = a,
    X = Xp + c1 exp(m1 x) + c2 exp(m2 x),   Xp = -a/(dS k^2),
    m1,2 = [-b +- sqrt(b^2 + 4 dS^2 k^2)] / (2 dS),

and c1, c2 from X(0) = X(Lx) = 0. dewi-lanjar evaluates this to machine
precision (analytic.exact_linear) and uses it to verify the numerics.
"""


@_topic("munk", "Munk (1950): lateral viscosity")
def _munk():
    return """
Steady, linear, viscous (dS = dI = 0):

    dM^3 lap(lap(psi)) - b psi_x = -F.

Boundary layer: with xi = x/dM, the balance is dM^3 psi_xxxx = b psi_x, so
psi_x ~ exp(m x) with

    dM^3 m^3 = b   =>   m = (b^(1/3)/dM) * {1, exp(+2 pi i/3), exp(-2 pi i/3)}.

The two complex roots have Re(m) = -b^(1/3)/(2 dM) < 0 and decay eastward:
an OSCILLATORY decay with wavenumber sqrt(3)/(2 dM). Two decaying roots let
the layer satisfy two wall conditions:

  no slip   (psi = psi_x = 0):
      psi = psi_I [1 - exp(-xi/2)(cos(sqrt3 xi/2) + sin(sqrt3 xi/2)/sqrt3)]
      v changes sign first at xi = 2 pi / sqrt3, i.e. x = 3.628 dM
  free slip (psi = psi_xx = 0):
      psi = psi_I [1 - exp(-xi/2)(cos(sqrt3 xi/2) - sin(sqrt3 xi/2)/sqrt3)]
      v changes sign first at xi = 4 pi / (3 sqrt3), i.e. x = 2.418 dM

The sign change of v is a weak offshore counter-current: a hallmark of the
Munk layer that the Stommel layer does not have.

With both drag and viscosity the local layer equation is
dM^3 m^3 - dS m - b = 0; dewi-lanjar uses its decaying roots for the
theoretical widths in every case (analytic.theory_widths). The exact
separable solution handles the full quartic in m and both wall conditions on
the western and eastern walls.
"""


@_topic("inertial", "Inertial boundary currents")
def _inertial():
    return """
When dI = (U/beta)^(1/2)/L is comparable with or larger than the frictional
widths, advection of relative vorticity, dI^2 J(psi, zeta), enters the
boundary-layer balance. In the inviscid limit the steady equation says that
potential vorticity q = dI^2 zeta + b y is constant along streamlines:

    J(psi, dI^2 zeta + b y) = 0.

A fluid column moving north in the western current must keep its q, so its
relative vorticity drops as it gains planetary vorticity. The current
cannot simply turn east where the linear current would: it overshoots
north and forms a recirculation cell, and the transport of the boundary
current exceeds the Sverdrup transport (transport_ratio > 1). Useful
ratios (reported by Params.regime()):

    Re_munk    = (dI/dM)^3,     Re_stommel = (dI/dS)^2.

For Re of order one or larger the boundary current is inertial; in a
double-gyre basin the merged currents separate as a free jet that may
become unstable and shed eddies.
"""


@_topic("fplane", "The f-plane control")
def _fplane():
    return """
Setting beta_factor = 0 removes the planetary vorticity gradient. The steady
drag-only problem becomes dS lap(psi) = F, whose solution for a wind that
does not depend on x is symmetric about x = Lx/2. There is no Sverdrup
interior and no western boundary current: intensification index = 1,
kinetic-energy centroid = 0.5. Running f_plane_control next to stommel_1948
isolates beta as the cause of the asymmetry; the sweep case beta_emergence
shows the transition continuously.
"""


@_topic("budgets", "Energy and vorticity budgets")
def _budgets():
    return """
Energy. Multiply the vorticity equation by -psi and integrate over the basin.
With psi = 0 on the walls,

    -integral psi zeta_t = d/dt [ (1/2) integral |grad psi|^2 ] = dKE/dt,
    integral psi J(psi, zeta) = 0,    integral psi psi_x = 0,

so

    dKE/dt = W - 2 dS KE - dM^3 integral zeta^2,     W = -integral psi F.

(The viscous boundary terms vanish for both free-slip and no-slip walls.)
In a steady state the wind work W equals the total dissipation. dewi-lanjar
evaluates each term with the discrete operators of the solver, so the
budget closes to round-off.

Vorticity. Integrate the equation over the basin. The beta term integrates
to b * integral (psi(Lx, y) - psi(0, y)) dy = 0: beta can only MOVE vorticity
around, not create it. In a steady state the wind input must therefore be
removed by friction, and the regional budgets show where:

    interior:       wind  ~ beta term            (Sverdrup)
    western strip:  beta  ~ drag + viscosity     (boundary layer)

On the discrete grid the basin-integrated beta term equals
b (psi_{N-1} - psi_1) / 2 per row, a boundary truncation error of O(h).
"""


@_topic("numerics", "Discretization and verification")
def _numerics():
    return """
Grid: uniform, node-centered; psi = 0 on the walls.

Space (second order):
  five-point Laplacian; centered d/dx; Arakawa (1966) Jacobian, which
  conserves discrete energy and enstrophy for advection and so removes
  nonlinear instability; wall vorticity zeta = 0 (free slip) or
  zeta_wall = 2 psi_1 / h^2 (no slip, Thom 1933).

Poisson problem lap(psi) = zeta: solved exactly in O(N log N) with the
type-I discrete sine transform, whose basis vectors are the eigenvectors of
the five-point Laplacian with Dirichlet data.

Steady linear problem: one sparse direct solve (SuperLU) of
[dM^3 (L L + B) - dS L - b Dx] psi = -F, with the SAME operators as the
time stepper, so a long linear integration converges to it to round-off.

Time: SSP-RK3 (Shu and Osher 1988). The step is the minimum of the
advective (sqrt3 / (dI^2 (|u|/dx + |v|/dy))), Rossby-wave
(sqrt3 / (b Ly / 2 pi)), viscous (2.51 / (dM^3 (4/dx^2 + 4/dy^2)) / 2) and
drag (2.51 / dS) limits times a safety factor.

Verification: against the exact separable solution, the error falls as h^2
(run `dewi-lanjar validate`). Resolution rule of thumb: at least three to
four grid points per boundary-layer decay scale; the diagnostics warn below
three.
"""


@_topic("diagnostics", "What the diagnostics mean")
def _diagnostics():
    return """
KE, Z                kinetic energy (1/2 int |u|^2) and enstrophy (1/2 int zeta^2).
energy budget        wind work, drag, viscous dissipation, advection, beta.
section_y            latitude of the gyre center (largest |psi|).
wbc_v_peak/x_peak    peak boundary-current velocity and its distance from the wall.
wbc_e_fold           distance from the wall at which |v| has decayed to |v_peak|/e.
wbc_zero_crossing    first sign change of v: edge of the counter-current.
theory_*             the same metrics for the boundary-layer theory profile.
wbc_transport        psi at the outer edge of the boundary current (where v first
                     changes sign): the current's volume transport.
transport_ratio      wbc_transport / Sverdrup psi_I at that edge; 1 when the
                     interior follows Sverdrup, < 1 when interior friction is not
                     negligible (large delta), > 1 when a recirculation adds transport.
intensification_index  peak |v| in the western half / peak |v| in the eastern half.
ke_centroid_x        zonal center of mass of kinetic energy, 0.5 if symmetric.
ke_entropy           Shannon entropy of the KE distribution / log N; 1 = uniform.
ke_zonal_entropy     same for the zonal marginal of KE: drops as energy
                     concentrates in the western boundary layer.
points_per_layer     grid points per boundary-layer decay scale (resolution).
grid_scale_fraction  share of the vorticity norm in two-grid-interval noise; above
                     about 0.3 the run is under-resolved.
vorticity_budget     wind, beta, drag, viscosity and advection integrated over
                     the basin, the gyre band, its western strip and interior.
"""


@_topic("floats", "Lagrangian floats")
def _floats():
    return """
Floats are passive particles moved by the velocity (u, v) = (-psi_y, psi_x):

    dx/dt = dI^2 u,    dy/dt = dI^2 v      (nondimensional).

The factor dI^2 = U/(beta L^2) converts the velocity unit U and the time
unit 1/(beta L) into a displacement in units of L. Linear runs drop the
advection of vorticity but water still moves, so floats there use a
kinematic delta_I set in [floats] delta_i. Integration: Heun's second-order
method with bilinear interpolation of the node velocities.
"""


def topic_names() -> list[str]:
    return list(TOPICS)


def explain(topic: str) -> str:
    key = topic.lower().replace("_", "-")
    aliases = {"nondim": "nondimensionalization", "scaling": "nondimensionalization",
               "f-plane": "fplane", "budget": "budgets", "numeric": "numerics",
               "diagnostic": "diagnostics", "float": "floats"}
    key = aliases.get(key, key)
    if key not in TOPICS:
        raise KeyError(f"unknown topic {topic!r}; choose from {', '.join(TOPICS)}")
    title, body = TOPICS[key]
    bar = "=" * len(title)
    return f"{title}\n{bar}\n{body.strip()}\n"
