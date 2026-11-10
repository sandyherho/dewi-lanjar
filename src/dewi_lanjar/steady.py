"""Steady linear solver: one sparse direct solve.

The steady linear problem

.. math::

    \\delta_M^3 \\nabla^2\\zeta - \\delta_S \\zeta - b\\,\\psi_x = -F ,\\qquad \\zeta = \\nabla^2\\psi

is discretized with exactly the operators of the time-dependent solver
(five-point Laplacian, centered ``d/dx``, wall vorticity from the boundary
condition). The interior unknowns ``psi`` obey ``M psi = -F`` with

``M = dM^3 (L L + B) - dS L - b Dx``

where ``L`` is the Dirichlet five-point Laplacian, ``Dx`` the centered
difference and ``B`` a diagonal matrix that carries the Thom wall vorticity
into the stencil of the first interior row (``B = 0`` for free slip). The
system is solved with SuperLU. Because the operators coincide, the time
integration of :mod:`dewi_lanjar.transient` relaxes to this solution to
round-off when ``delta_i = 0``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .forcing import WindForcing
from .grid import Grid
from .operators import laplacian, wall_vorticity
from .params import Params


def _second_difference(n: int, h: float) -> sp.csr_matrix:
    m = n - 1
    return sp.diags([np.ones(m - 1), -2 * np.ones(m), np.ones(m - 1)], [-1, 0, 1]) / h**2


def _first_difference(n: int, h: float) -> sp.csr_matrix:
    m = n - 1
    return sp.diags([-np.ones(m - 1), np.ones(m - 1)], [-1, 1]) / (2 * h)


def assemble_operators(grid: Grid, bc: str = "no-slip", bc_ns: str | None = None) -> dict:
    """Sparse interior operators (row-major, index ``i * (ny-1) + j``)."""
    Ix = sp.identity(grid.nx - 1)
    Iy = sp.identity(grid.ny - 1)
    L = sp.kron(_second_difference(grid.nx, grid.dx), Iy) + sp.kron(Ix, _second_difference(grid.ny, grid.dy))
    Dx = sp.kron(_first_difference(grid.nx, grid.dx), Iy)
    bc_ns = bc if bc_ns is None else bc_ns
    bx = np.zeros(grid.interior_shape)
    if bc == "no-slip":
        bx[0, :] += 2.0 / grid.dx**4
        bx[-1, :] += 2.0 / grid.dx**4
    if bc_ns == "no-slip":
        bx[:, 0] += 2.0 / grid.dy**4
        bx[:, -1] += 2.0 / grid.dy**4
    B = sp.diags(bx.ravel())
    return {"L": L.tocsr(), "Dx": Dx.tocsr(), "B": B.tocsr()}


@dataclass
class SteadySolution:
    """Result of :func:`solve_steady`."""

    grid: Grid
    params: Params
    forcing: WindForcing
    psi: np.ndarray
    zeta: np.ndarray
    residual: float
    info: dict = field(default_factory=dict)


def solve_steady(grid: Grid, p: Params, forcing: WindForcing) -> SteadySolution:
    """Solve the steady *linear* problem on ``grid`` (``delta_i`` is ignored).

    Raises
    ------
    ValueError
        If neither drag nor viscosity is present (no steady solution exists).
    """
    if not p.has_friction:
        raise ValueError("a steady state needs delta_s > 0 or delta_m > 0")
    ops = assemble_operators(grid, p.bc, p.bc_zonal_walls)
    L, Dx, B = ops["L"], ops["Dx"], ops["B"]
    M = -p.delta_s * L - p.beta_factor * Dx
    if p.dm3 > 0:
        M = M + p.dm3 * (L @ L + B)
    F = forcing.curl_steady(grid.X, grid.Y)[1:-1, 1:-1].ravel()
    psi_int = spla.spsolve(M.tocsc(), -F)
    res = float(np.max(np.abs(M @ psi_int + F)) / max(np.max(np.abs(F)), 1e-300))
    psi = grid.pad(psi_int.reshape(grid.interior_shape))
    zeta = grid.pad(laplacian(psi, grid.dx, grid.dy))
    wall_vorticity(psi, zeta, grid.dx, grid.dy, p.bc, p.bc_zonal_walls)
    return SteadySolution(grid=grid, params=p, forcing=forcing, psi=psi, zeta=zeta,
                          residual=res, info={"n_unknowns": psi_int.size, "nnz": M.nnz})
