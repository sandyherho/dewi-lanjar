"""Uniform node-centered grid on the rectangle [0, Lx] x [0, Ly].

Arrays are indexed ``[i, j]`` with ``i`` along x (west to east) and ``j``
along y (south to north). Node ``(0, j)`` lies on the western wall. The
streamfunction vanishes on all four walls, so only the
``(nx - 1) x (ny - 1)`` interior nodes carry unknowns.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz


@dataclass(frozen=True)
class Grid:
    """Uniform grid with ``nx`` x ``ny`` cells (``(nx+1) x (ny+1)`` nodes)."""

    nx: int = 128
    ny: int = 128
    Lx: float = 1.0
    Ly: float = 1.0

    def __post_init__(self) -> None:
        if int(self.nx) != self.nx or int(self.ny) != self.ny:
            raise ValueError("nx and ny must be integers")
        if self.nx < 4 or self.ny < 4:
            raise ValueError("nx and ny must be at least 4")
        if not (self.Lx > 0 and self.Ly > 0):
            raise ValueError("Lx and Ly must be positive")

    @property
    def dx(self) -> float:
        return self.Lx / self.nx

    @property
    def dy(self) -> float:
        return self.Ly / self.ny

    @property
    def h(self) -> float:
        """Largest grid spacing."""
        return max(self.dx, self.dy)

    @cached_property
    def x(self) -> np.ndarray:
        return np.linspace(0.0, self.Lx, self.nx + 1)

    @cached_property
    def y(self) -> np.ndarray:
        return np.linspace(0.0, self.Ly, self.ny + 1)

    @cached_property
    def X(self) -> np.ndarray:
        return np.meshgrid(self.x, self.y, indexing="ij")[0]

    @cached_property
    def Y(self) -> np.ndarray:
        return np.meshgrid(self.x, self.y, indexing="ij")[1]

    @property
    def shape(self) -> tuple[int, int]:
        return (self.nx + 1, self.ny + 1)

    @property
    def interior_shape(self) -> tuple[int, int]:
        return (self.nx - 1, self.ny - 1)

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy

    def integrate(self, f: np.ndarray) -> float:
        """Trapezoidal integral of a full-grid field over the basin."""
        return float(_trapz(_trapz(f, self.y, axis=1), self.x))

    def pad(self, interior: np.ndarray, fill: float = 0.0) -> np.ndarray:
        """Embed an interior array in a full array with constant walls."""
        out = np.full(self.shape, fill, dtype=interior.dtype)
        out[1:-1, 1:-1] = interior
        return out
