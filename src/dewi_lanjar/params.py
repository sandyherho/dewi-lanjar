"""Model parameters, nondimensionalization and boundary-layer scales.

dewi-lanjar solves the barotropic vorticity equation on a beta plane in a
closed rectangular basin. Every run works in the nondimensional form

.. math::

    \\partial_t \\zeta + \\delta_I^2 J(\\psi, \\zeta) + b\\,\\partial_x \\psi
    = F(x, y, t) - \\delta_S \\zeta + \\delta_M^3 \\nabla^2 \\zeta ,
    \\qquad \\zeta = \\nabla^2 \\psi ,

in which lengths are scaled by the zonal basin width ``L``, velocities by
the Sverdrup velocity ``U = tau0 / (rho0 H beta L)`` and time by
``1 / (beta L)``. The three nondimensional numbers are the boundary-layer
widths divided by ``L``:

* ``delta_S = r / (beta L)``                (Stommel, bottom drag)
* ``delta_M = (A_H / beta)**(1/3) / L``     (Munk, lateral viscosity)
* ``delta_I = (U / beta)**(1/2) / L``       (inertial, advection)

``b`` (``beta_factor``) multiplies the planetary vorticity gradient. It
equals 1 for a beta plane and 0 for an f plane, the control experiment in
which western intensification disappears.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Optional

VALID_BC = ("free-slip", "no-slip")

#: Earth's rotation rate [rad/s].
OMEGA_EARTH = 7.2921159e-5
#: Mean Earth radius [m].
R_EARTH = 6.371e6


def beta_at_latitude(lat_deg: float) -> float:
    """Planetary vorticity gradient ``beta = 2 Omega cos(phi) / a`` [1/(m s)]."""
    return 2.0 * OMEGA_EARTH * math.cos(math.radians(lat_deg)) / R_EARTH


@dataclass(frozen=True)
class Params:
    """Nondimensional parameters of the barotropic vorticity model.

    Parameters
    ----------
    delta_s, delta_m, delta_i:
        Stommel, Munk and inertial boundary-layer widths divided by the
        basin width. Each must be non-negative.
    beta_factor:
        Multiplier ``b`` of the beta term (1 = beta plane, 0 = f plane).
    bc:
        Lateral boundary condition on the western and eastern walls,
        ``"free-slip"`` (zero vorticity) or ``"no-slip"`` (zero tangential
        velocity). It matters only when ``delta_m > 0``.
    bc_ns:
        Condition on the southern and northern walls. ``None`` (default)
        uses ``bc``. The exact separable solution of
        :mod:`dewi_lanjar.analytic` requires ``"free-slip"`` here.
    Lx, Ly:
        Basin size in units of ``L``. ``Lx`` is normally 1.
    """

    delta_s: float = 0.0
    delta_m: float = 0.0
    delta_i: float = 0.0
    beta_factor: float = 1.0
    bc: str = "free-slip"
    Lx: float = 1.0
    Ly: float = 1.0
    bc_ns: Optional[str] = None

    def __post_init__(self) -> None:
        for name in ("delta_s", "delta_m", "delta_i"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v >= 0):
                raise ValueError(f"{name} must be a finite number >= 0, got {v!r}")
        if not math.isfinite(self.beta_factor):
            raise ValueError("beta_factor must be finite")
        if self.bc not in VALID_BC:
            raise ValueError(f"bc must be one of {VALID_BC}, got {self.bc!r}")
        if self.bc_ns is not None and self.bc_ns not in VALID_BC:
            raise ValueError(f"bc_ns must be None or one of {VALID_BC}, got {self.bc_ns!r}")
        if not (self.Lx > 0 and self.Ly > 0):
            raise ValueError("Lx and Ly must be positive")

    # ------------------------------------------------------------------
    @property
    def bc_zonal_walls(self) -> str:
        """Effective condition on the southern and northern walls."""
        return self.bc if self.bc_ns is None else self.bc_ns

    @property
    def dm3(self) -> float:
        """``delta_M**3``, the coefficient of the viscous term."""
        return self.delta_m ** 3

    @property
    def di2(self) -> float:
        """``delta_I**2``, the coefficient of the advective term."""
        return self.delta_i ** 2

    @property
    def is_linear(self) -> bool:
        return self.delta_i == 0.0

    @property
    def has_friction(self) -> bool:
        return self.delta_s > 0.0 or self.delta_m > 0.0

    def with_(self, **kw) -> "Params":
        """Return a copy with some fields replaced."""
        return replace(self, **kw)

    # ------------------------------------------------------------------
    def regime(self) -> dict:
        """Dimensionless ratios that classify the dynamical regime.

        ``inertial_stommel = delta_I / delta_S`` and
        ``inertial_munk = delta_I / delta_M`` measure the strength of
        advection against friction inside the western boundary layer: a
        boundary-layer Reynolds number is ``(delta_I/delta_M)**3`` for
        viscous closure and ``(delta_I/delta_S)**2`` for drag closure.
        Values well above one mean that the boundary current is inertial
        (overshoot, recirculation, possibly eddies).
        """
        out = {
            "delta_s": self.delta_s,
            "delta_m": self.delta_m,
            "delta_i": self.delta_i,
            "beta_factor": self.beta_factor,
        }
        out["inertial_stommel"] = (
            self.delta_i / self.delta_s if self.delta_s > 0 else math.inf
        )
        out["inertial_munk"] = (
            self.delta_i / self.delta_m if self.delta_m > 0 else math.inf
        )
        out["Re_munk"] = (
            (self.delta_i / self.delta_m) ** 3 if self.delta_m > 0 else math.inf
        )
        out["Re_stommel"] = (
            (self.delta_i / self.delta_s) ** 2 if self.delta_s > 0 else math.inf
        )
        if self.beta_factor == 0:
            label = "f-plane (no planetary vorticity gradient)"
        elif self.delta_i == 0:
            if self.delta_m == 0:
                label = "linear Stommel (drag-closed)"
            elif self.delta_s == 0:
                label = "linear Munk (viscous)"
            else:
                label = "linear Stommel-Munk"
        else:
            r = max(
                self.delta_i / self.delta_m if self.delta_m > 0 else 0.0,
                self.delta_i / self.delta_s if self.delta_s > 0 else 0.0,
            )
            label = (
                "weakly nonlinear" if r < 1.0 else "inertial (nonlinear boundary current)"
            )
        out["label"] = label
        return out


@dataclass(frozen=True)
class DimensionalScales:
    """Dimensional parameters and the scales that map them to :class:`Params`.

    All inputs are SI. ``L`` is the zonal basin width.

    Parameters
    ----------
    L:
        Zonal basin width [m].
    beta:
        Planetary vorticity gradient [1/(m s)].
    r:
        Linear bottom-drag coefficient [1/s].
    A_h:
        Lateral eddy viscosity [m^2/s].
    tau0:
        Wind-stress amplitude [N/m^2].
    rho0:
        Reference density [kg/m^3].
    H:
        Layer depth [m].
    """

    L: float = 4.0e6
    beta: float = 2.0e-11
    r: float = 1.0e-6
    A_h: float = 0.0
    tau0: float = 0.1
    rho0: float = 1025.0
    H: float = 1000.0
    _extra: dict = field(default_factory=dict, compare=False, repr=False)

    def __post_init__(self) -> None:
        for name in ("L", "beta", "tau0", "rho0", "H"):
            if not getattr(self, name) > 0:
                raise ValueError(f"{name} must be positive")
        for name in ("r", "A_h"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0")

    # scales -----------------------------------------------------------
    @property
    def U(self) -> float:
        """Sverdrup velocity scale ``tau0 / (rho0 H beta L)`` [m/s]."""
        return self.tau0 / (self.rho0 * self.H * self.beta * self.L)

    @property
    def time_scale(self) -> float:
        """Time scale ``1 / (beta L)`` [s]."""
        return 1.0 / (self.beta * self.L)

    @property
    def transport_scale(self) -> float:
        """Volume transport per unit nondimensional streamfunction [m^3/s].

        ``H U L = tau0 / (rho0 beta)``; independent of ``H`` and ``L``.
        """
        return self.tau0 / (self.rho0 * self.beta)

    @property
    def delta_s(self) -> float:
        return self.r / (self.beta * self.L)

    @property
    def delta_m(self) -> float:
        return (self.A_h / self.beta) ** (1.0 / 3.0) / self.L

    @property
    def delta_i(self) -> float:
        return math.sqrt(self.U / self.beta) / self.L

    def to_params(self, Ly_over_L: float = 1.0, bc: str = "free-slip",
                  beta_factor: float = 1.0, inertia: bool = True) -> Params:
        """Build the nondimensional :class:`Params`."""
        return Params(
            delta_s=self.delta_s,
            delta_m=self.delta_m,
            delta_i=self.delta_i if inertia else 0.0,
            beta_factor=beta_factor,
            bc=bc,
            Lx=1.0,
            Ly=Ly_over_L,
        )

    # conversions -------------------------------------------------------
    def sverdrups(self, psi_nd: float) -> float:
        """Convert a nondimensional streamfunction to sverdrups (1e6 m^3/s)."""
        return psi_nd * self.transport_scale / 1.0e6

    def km(self, length_nd: float) -> float:
        return length_nd * self.L / 1.0e3

    def m_per_s(self, vel_nd: float) -> float:
        return vel_nd * self.U

    def days(self, t_nd: float) -> float:
        return t_nd * self.time_scale / 86400.0

    def explain(self) -> str:
        """Return a worked derivation of every scale with numbers substituted."""
        b, L = self.beta, self.L
        lines = [
            "Dimensional inputs",
            f"  L     = {L:.4g} m        (zonal basin width)",
            f"  beta  = {b:.4g} 1/(m s)  (planetary vorticity gradient)",
            f"  r     = {self.r:.4g} 1/s      (linear bottom drag)",
            f"  A_h   = {self.A_h:.4g} m^2/s    (lateral eddy viscosity)",
            f"  tau0  = {self.tau0:.4g} N/m^2    (wind-stress amplitude)",
            f"  rho0  = {self.rho0:.4g} kg/m^3",
            f"  H     = {self.H:.4g} m",
            "",
            "Sverdrup velocity scale (wind-curl input balanced by beta*v):",
            f"  U = tau0/(rho0 H beta L) = {self.tau0:.4g}/({self.rho0:.4g}*{self.H:.4g}"
            f"*{b:.4g}*{L:.4g}) = {self.U:.4g} m/s",
            "Time scale (inverse Rossby-wave frequency at the basin scale):",
            f"  T = 1/(beta L) = {self.time_scale:.4g} s = {self.time_scale / 86400:.4g} days",
            "Transport scale (H U L, independent of H):",
            f"  tau0/(rho0 beta) = {self.transport_scale:.4g} m^3/s"
            f" = {self.transport_scale / 1e6:.4g} Sv",
            "",
            "Boundary-layer widths (balance of beta*v against each term):",
            f"  Stommel  d_S = r/beta            = {self.r / b:.4g} m"
            f"  -> delta_S = {self.delta_s:.4g}",
            f"  Munk     d_M = (A_h/beta)^(1/3)  = {(self.A_h / b) ** (1 / 3):.4g} m"
            f"  -> delta_M = {self.delta_m:.4g}",
            f"  inertial d_I = (U/beta)^(1/2)    = {math.sqrt(self.U / b):.4g} m"
            f"  -> delta_I = {self.delta_i:.4g}",
        ]
        return "\n".join(lines)


def boundary_layer_roots(p: Params) -> list[complex]:
    """Roots of the local western-boundary-layer equation.

    Near the western wall the along-wall scale is long compared with the
    cross-wall scale, so the steady linear equation reduces to the ODE
    ``delta_M^3 X''' - delta_S X' - b X = 0`` for ``X = d(psi)/dx``,
    equivalently the cubic ``delta_M^3 m^3 - delta_S m - b = 0`` for
    solutions ``exp(m x)``. Roots with negative real part decay eastward
    from the western wall and form the western boundary current.
    """
    import numpy as np

    if p.dm3 > 0:
        roots = np.roots([p.dm3, 0.0, -p.delta_s, -p.beta_factor])
    elif p.delta_s > 0:
        roots = np.array([-p.beta_factor / p.delta_s])
    else:
        raise ValueError("a boundary layer needs delta_s > 0 or delta_m > 0")
    return [complex(z) for z in roots]
