"""Lagrangian floats (passive particles) advected by the model velocity.

In nondimensional variables a particle moves as

.. math::

    \\frac{d\\mathbf{x}}{dt} = \\delta_I^2\\,\\mathbf{u}(\\mathbf{x}, t),
    \\qquad \\mathbf{u} = (-\\psi_y, \\psi_x),

because the time unit is ``1/(beta L)`` and the velocity unit is the
Sverdrup velocity ``U``. The factor ``delta_I^2 = U / (beta L^2)`` is the
same factor that multiplies the advection of vorticity. A *linear* model
drops vorticity advection but the water still moves, so floats take their
own kinematic ``delta_I`` (``kinematic_delta_i``) that need not match the
dynamical one.

Velocities are interpolated bilinearly from the grid nodes; the time step
is Heun's (explicit trapezoidal) second-order scheme. Floats that reach a
wall are reflected back inside (the wall-normal velocity vanishes there, so
this only matters at the level of round-off and time-step error).
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates

from .grid import Grid
from .operators import velocity


class Floats:
    """A cloud of passive particles."""

    def __init__(self, grid: Grid, n: int, seed: int = 0, speed_factor: float = 1.0,
                 margin: float = 0.01) -> None:
        rng = np.random.default_rng(seed)
        self.grid = grid
        self.speed_factor = float(speed_factor)
        mx, my = margin * grid.Lx, margin * grid.Ly
        self.pos = np.column_stack([
            rng.uniform(mx, grid.Lx - mx, n),
            rng.uniform(my, grid.Ly - my, n),
        ])

    def _interp(self, field: np.ndarray, pos: np.ndarray) -> np.ndarray:
        g = self.grid
        coords = np.vstack([pos[:, 0] / g.dx, pos[:, 1] / g.dy])
        return map_coordinates(field, coords, order=1, mode="nearest")

    def velocity_at(self, u: np.ndarray, v: np.ndarray, pos: np.ndarray) -> np.ndarray:
        return self.speed_factor * np.column_stack([self._interp(u, pos), self._interp(v, pos)])

    def _reflect(self, pos: np.ndarray) -> np.ndarray:
        g = self.grid
        x, y = pos[:, 0], pos[:, 1]
        x = np.where(x < 0, -x, x)
        x = np.where(x > g.Lx, 2 * g.Lx - x, x)
        y = np.where(y < 0, -y, y)
        y = np.where(y > g.Ly, 2 * g.Ly - y, y)
        return np.column_stack([np.clip(x, 0, g.Lx), np.clip(y, 0, g.Ly)])

    def step(self, psi_old: np.ndarray, psi_new: np.ndarray, dt: float) -> None:
        """Advance by ``dt`` with the flow changing from ``psi_old`` to ``psi_new``."""
        g = self.grid
        u0, v0 = velocity(psi_old, g.dx, g.dy)
        if psi_new is psi_old:
            u1, v1 = u0, v0
        else:
            u1, v1 = velocity(psi_new, g.dx, g.dy)
        k1 = self.velocity_at(u0, v0, self.pos)
        pred = self._reflect(self.pos + dt * k1)
        k2 = self.velocity_at(u1, v1, pred)
        self.pos = self._reflect(self.pos + 0.5 * dt * (k1 + k2))

    def advect_steady(self, psi: np.ndarray, t_total: float, n_frames: int,
                      substeps: int = 8) -> np.ndarray:
        """Advect through a steady flow; return positions, shape ``(n_frames, n, 2)``."""
        g = self.grid
        u, v = velocity(psi, g.dx, g.dy)
        dt_frame = t_total / max(n_frames - 1, 1)
        smax = self.speed_factor * float(np.max(np.hypot(u, v)))
        # keep each substep below half a grid cell of travel
        substeps = max(substeps, int(np.ceil(smax * dt_frame / (0.5 * min(g.dx, g.dy)))))
        dt = dt_frame / substeps
        out = np.empty((n_frames,) + self.pos.shape)
        out[0] = self.pos
        for f in range(1, n_frames):
            for _ in range(substeps):
                k1 = self.velocity_at(u, v, self.pos)
                pred = self._reflect(self.pos + dt * k1)
                k2 = self.velocity_at(u, v, pred)
                self.pos = self._reflect(self.pos + 0.5 * dt * (k1 + k2))
            out[f] = self.pos
        return out
