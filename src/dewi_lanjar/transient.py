"""Time-dependent nonlinear solver.

Integrates

.. math::

    \\partial_t\\zeta = -\\delta_I^2 J(\\psi,\\zeta) - b\\,\\psi_x + F(x,y,t)
    - \\delta_S\\zeta + \\delta_M^3\\nabla^2\\zeta ,\\qquad \\nabla^2\\psi = \\zeta

on the interior nodes with

* the Arakawa (1966) Jacobian (discrete energy conservation);
* the exact discrete-sine-transform Poisson inversion;
* wall vorticity from the lateral boundary condition, refreshed at every
  stage;
* the three-stage, third-order strong-stability-preserving Runge-Kutta
  scheme of Shu and Osher (1988);
* an adaptive time step limited by advection, Rossby-wave, viscous and
  drag stability bounds, trimmed so that output times are hit exactly.
"""

from __future__ import annotations

import math
import time as _time
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from . import diagnostics as diag
from .floats import Floats
from .forcing import WindForcing
from .grid import Grid
from .operators import PoissonSolver, arakawa_jacobian, ddx, laplacian, velocity, wall_vorticity
from .params import Params

# stability radii of SSP-RK3 along the imaginary and negative real axes
_RK3_IMAG = math.sqrt(3.0)
_RK3_REAL = 2.5127


@dataclass
class RunConfig:
    """Controls a time integration.

    Parameters
    ----------
    t_end:
        Final nondimensional time (unit ``1/(beta L)``).
    n_frames:
        Number of stored snapshots, evenly spaced on ``[t_record, t_end]``.
    t_record:
        Time of the first stored snapshot (skip the spin-up if desired).
    cfl:
        Safety factor (0 < cfl <= 1) applied to every stability bound.
    dt_max:
        Optional upper limit on the time step.
    n_floats:
        Number of Lagrangian floats (0 for none).
    float_seed:
        Random seed for the initial float positions.
    float_delta_i:
        Kinematic inertial width used to move floats. ``None`` uses the
        dynamical ``delta_i`` of the parameters, or 0.03 for a linear run.
    diag_per_frame:
        Diagnostics are recorded this many times per frame interval.
    init_noise:
        Amplitude (relative to the Sverdrup vorticity scale, which is 1) of a
        smooth random initial vorticity field. A small value breaks exact
        symmetries of the forcing that the discretization would otherwise
        preserve forever (for example in a symmetric double gyre).
    init_seed:
        Seed of the initial noise.
    store_dtype:
        Precision of stored snapshots (``float32`` halves memory).
    """

    t_end: float = 100.0
    n_frames: int = 120
    t_record: float = 0.0
    cfl: float = 0.6
    dt_max: Optional[float] = None
    n_floats: int = 0
    float_seed: int = 0
    float_delta_i: Optional[float] = None
    diag_per_frame: int = 2
    init_noise: float = 0.0
    init_seed: int = 0
    store_dtype: str = "float32"

    def __post_init__(self) -> None:
        if not self.t_end > 0:
            raise ValueError("t_end must be positive")
        if not 0 <= self.t_record < self.t_end:
            raise ValueError("t_record must lie in [0, t_end)")
        if self.n_frames < 2:
            raise ValueError("n_frames must be at least 2")
        if not 0 < self.cfl <= 1:
            raise ValueError("cfl must lie in (0, 1]")
        if self.n_floats < 0:
            raise ValueError("n_floats must be >= 0")


@dataclass
class TransientResult:
    """Output of :meth:`TransientSolver.run`."""

    grid: Grid
    params: Params
    forcing: WindForcing
    config: RunConfig
    times: np.ndarray
    psi: np.ndarray            # (n_frames, nx+1, ny+1)
    zeta: np.ndarray           # (n_frames, nx+1, ny+1)
    floats: Optional[np.ndarray]   # (n_frames, n_floats, 2) or None
    series: dict               # time series of scalar diagnostics
    final: dict                # full diagnostics of the last state
    stats: dict = field(default_factory=dict)

    @property
    def psi_final(self) -> np.ndarray:
        return self.psi[-1].astype(float)


class TransientSolver:
    """Time integrator for the barotropic vorticity equation."""

    def __init__(self, grid: Grid, params: Params, forcing: WindForcing) -> None:
        if abs(grid.Lx - params.Lx) > 1e-12 or abs(grid.Ly - params.Ly) > 1e-12:
            raise ValueError("grid and params describe different basins")
        self.grid = grid
        self.p = params
        self.forcing = forcing
        self.poisson = PoissonSolver(grid.nx, grid.ny, grid.dx, grid.dy)
        self._F0 = forcing.curl_steady(grid.X, grid.Y)[1:-1, 1:-1]

    # ------------------------------------------------------------------
    def fields(self, zeta_int: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Full ``psi`` and ``zeta`` (walls included) from interior vorticity."""
        g = self.grid
        psi = g.pad(self.poisson.solve(zeta_int))
        zeta = g.pad(zeta_int)
        wall_vorticity(psi, zeta, g.dx, g.dy, self.p.bc, self.p.bc_zonal_walls)
        return psi, zeta

    def rhs(self, zeta_int: np.ndarray, t: float) -> np.ndarray:
        """Right-hand side ``d zeta / dt`` at the interior nodes."""
        g, p = self.grid, self.p
        psi, zeta = self.fields(zeta_int)
        out = self._F0 * self.forcing.modulation(t)
        if p.beta_factor != 0:
            out = out - p.beta_factor * ddx(psi, g.dx)
        if p.delta_s > 0:
            out = out - p.delta_s * zeta_int
        if p.dm3 > 0:
            out = out + p.dm3 * laplacian(zeta, g.dx, g.dy)
        if p.di2 > 0:
            out = out - p.di2 * arakawa_jacobian(psi, zeta, g.dx, g.dy)
        return out

    def stable_dt(self, psi: np.ndarray, cfl: float) -> float:
        """Largest stable time step for the current state.

        * advection: ``dt <= sqrt(3) / (dI^2 (|u|/dx + |v|/dy))``;
        * Rossby waves: ``omega <= b k_x / k^2 <= b Ly / (2 pi)``;
        * viscosity: ``dt <= 2.51 / (dM^3 (4/dx^2 + 4/dy^2))``;
        * drag: ``dt <= 2.51 / dS``.
        """
        g, p = self.grid, self.p
        bounds = []
        if p.di2 > 0:
            u, v = velocity(psi, g.dx, g.dy)
            rate = p.di2 * float(np.max(np.abs(u) / g.dx + np.abs(v) / g.dy))
            if rate > 0:
                bounds.append(_RK3_IMAG / rate)
        if p.beta_factor != 0:
            bounds.append(_RK3_IMAG / (abs(p.beta_factor) * g.Ly / (2 * math.pi)))
        if p.dm3 > 0:
            # no-slip walls add stiffness; the factor 0.5 covers the Thom closure
            bounds.append(0.5 * _RK3_REAL / (p.dm3 * (4 / g.dx**2 + 4 / g.dy**2)))
        if p.delta_s > 0:
            bounds.append(_RK3_REAL / p.delta_s)
        dt = cfl * min(bounds) if bounds else cfl * 1.0
        return dt

    def step(self, zeta_int: np.ndarray, t: float, dt: float) -> np.ndarray:
        """One SSP-RK3 step (Shu and Osher 1988)."""
        z1 = zeta_int + dt * self.rhs(zeta_int, t)
        z2 = 0.75 * zeta_int + 0.25 * (z1 + dt * self.rhs(z1, t + dt))
        return zeta_int / 3.0 + 2.0 / 3.0 * (z2 + dt * self.rhs(z2, t + 0.5 * dt))

    # ------------------------------------------------------------------
    def run(self, config: RunConfig, zeta0: Optional[np.ndarray] = None,
            progress: Optional[Callable[[float, float], None]] = None) -> TransientResult:
        """Integrate from rest (or from ``zeta0``, full or interior) to ``config.t_end``."""
        g, p = self.grid, self.p
        if zeta0 is None:
            z = np.zeros(g.interior_shape)
            if config.init_noise > 0:
                from scipy.ndimage import gaussian_filter
                rng = np.random.default_rng(config.init_seed)
                noise = gaussian_filter(rng.standard_normal(g.interior_shape), 3.0)
                z = config.init_noise * noise / np.max(np.abs(noise))
        else:
            z = np.array(zeta0, float)
            if z.shape == g.shape:
                z = z[1:-1, 1:-1].copy()
            if z.shape != g.interior_shape:
                raise ValueError("zeta0 has the wrong shape")

        frame_times = np.linspace(config.t_record, config.t_end, config.n_frames)
        n_diag = (config.n_frames - 1) * config.diag_per_frame + 1
        diag_times = np.linspace(config.t_record, config.t_end, n_diag)
        if config.t_record > 0:
            pre = np.linspace(0.0, config.t_record, max(2, config.diag_per_frame * 10))[:-1]
            diag_times = np.concatenate([pre, diag_times])
        events = np.unique(np.concatenate([frame_times, diag_times]))
        events = events[events > 0] if events[0] == 0 else events

        dtype = np.dtype(config.store_dtype)
        psi_store = np.empty((config.n_frames,) + g.shape, dtype)
        zeta_store = np.empty((config.n_frames,) + g.shape, dtype)
        floats = None
        float_store = None
        if config.n_floats > 0:
            kdi = config.float_delta_i
            if kdi is None:
                kdi = p.delta_i if p.delta_i > 0 else 0.03
            floats = Floats(g, config.n_floats, seed=config.float_seed, speed_factor=kdi**2)
            float_store = np.empty((config.n_frames, config.n_floats, 2))

        series: dict = {k: [] for k in (
            "t", "KE", "enstrophy", "psi_max", "psi_min", "wbc_transport", "transport_ratio",
            "intensification_index", "ke_entropy", "ke_zonal_entropy", "ke_centroid_x",
            "energy_wind", "energy_drag", "energy_viscous", "dKE_dt", "max_speed")}

        def record_series(psi, zeta, t):
            d = diag.compute(psi, zeta, g, p, self.forcing, t=t, full=False)
            for k in series:
                series[k].append(d[k])

        def record_frame(idx, psi, zeta):
            psi_store[idx] = psi
            zeta_store[idx] = zeta
            if floats is not None:
                float_store[idx] = floats.pos

        t = 0.0
        psi, zeta = self.fields(z)
        fi = 0
        di_set = set(np.round(diag_times, 12))
        if frame_times[0] == 0.0:
            record_frame(0, psi, zeta)
            fi = 1
        if 0.0 in di_set:
            record_series(psi, zeta, 0.0)

        nsteps = 0
        wall0 = _time.perf_counter()
        dt_min, dt_maxseen = math.inf, 0.0
        for target in events:
            while t < target - 1e-12:
                dt = self.stable_dt(psi, config.cfl)
                if config.dt_max is not None:
                    dt = min(dt, config.dt_max)
                dt = min(dt, target - t)
                psi_old = psi
                z = self.step(z, t, dt)
                t += dt
                nsteps += 1
                dt_min, dt_maxseen = min(dt_min, dt), max(dt_maxseen, dt)
                psi, zeta = self.fields(z)
                if floats is not None:
                    floats.step(psi_old, psi, dt)
                if not np.isfinite(z).all():
                    raise FloatingPointError(
                        f"solution blew up at t = {t:.4g}; refine the grid or reduce cfl"
                    )
            t = float(target)
            if fi < config.n_frames and abs(t - frame_times[fi]) < 1e-9:
                record_frame(fi, psi, zeta)
                fi += 1
            if round(t, 12) in di_set:
                record_series(psi, zeta, t)
            if progress is not None:
                progress(t, config.t_end)

        final = diag.compute(psi, zeta, g, p, self.forcing, t=t, full=True)
        stats = {
            "steps": nsteps,
            "dt_min": dt_min,
            "dt_max": dt_maxseen,
            "wall_seconds": _time.perf_counter() - wall0,
        }
        return TransientResult(
            grid=g, params=p, forcing=self.forcing, config=config,
            times=frame_times, psi=psi_store, zeta=zeta_store, floats=float_store,
            series={k: np.asarray(v) for k, v in series.items()}, final=final, stats=stats,
        )
