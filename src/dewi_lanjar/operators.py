"""Second-order finite-difference operators on the node-centered grid.

All operators act on *full* arrays of shape ``(nx+1, ny+1)`` and return
*interior* arrays of shape ``(nx-1, ny-1)``.

* :func:`laplacian` is the five-point Laplacian.
* :func:`ddx`, :func:`ddy` are centered first differences.
* :func:`arakawa_jacobian` is the Arakawa (1966) nine-point Jacobian. With
  ``psi`` constant on the boundary it conserves the discrete kinetic
  energy ``-1/2 sum(psi zeta)`` exactly, so advection can neither create
  nor destroy energy, which prevents nonlinear instability.
* :class:`PoissonSolver` inverts the five-point Laplacian with homogeneous
  Dirichlet data exactly (to round-off) through the type-I discrete sine
  transform, whose basis vectors are the eigenvectors of that matrix.
* :func:`wall_vorticity` sets the vorticity on the walls from the lateral
  boundary condition (zero for free slip, the Thom formula for no slip).
"""

from __future__ import annotations

import numpy as np
from scipy import fft as sfft


def laplacian(f: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Five-point Laplacian at interior nodes."""
    return (
        (f[2:, 1:-1] - 2.0 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / dx**2
        + (f[1:-1, 2:] - 2.0 * f[1:-1, 1:-1] + f[1:-1, :-2]) / dy**2
    )


def ddx(f: np.ndarray, dx: float) -> np.ndarray:
    """Centered x derivative at interior nodes."""
    return (f[2:, 1:-1] - f[:-2, 1:-1]) / (2.0 * dx)


def ddy(f: np.ndarray, dy: float) -> np.ndarray:
    """Centered y derivative at interior nodes."""
    return (f[1:-1, 2:] - f[1:-1, :-2]) / (2.0 * dy)


def arakawa_jacobian(psi: np.ndarray, zeta: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Arakawa (1966) Jacobian ``J(psi, zeta) = psi_x zeta_y - psi_y zeta_x``.

    The average of three second-order forms (J++, J+x, Jx+). The average
    inherits the discrete analogues of
    ``int psi J = 0`` (energy) and ``int zeta J = 0`` (enstrophy).
    """
    p, z = psi, zeta
    # shifted views
    pE, pW, pN, pS = p[2:, 1:-1], p[:-2, 1:-1], p[1:-1, 2:], p[1:-1, :-2]
    pNE, pNW, pSE, pSW = p[2:, 2:], p[:-2, 2:], p[2:, :-2], p[:-2, :-2]
    zE, zW, zN, zS = z[2:, 1:-1], z[:-2, 1:-1], z[1:-1, 2:], z[1:-1, :-2]
    zNE, zNW, zSE, zSW = z[2:, 2:], z[:-2, 2:], z[2:, :-2], z[:-2, :-2]

    jpp = (pE - pW) * (zN - zS) - (pN - pS) * (zE - zW)
    jpx = (
        pE * (zNE - zSE) - pW * (zNW - zSW)
        - pN * (zNE - zNW) + pS * (zSE - zSW)
    )
    jxp = (
        zN * (pNE - pNW) - zS * (pSE - pSW)
        - zE * (pNE - pSE) + zW * (pNW - pSW)
    )
    return (jpp + jpx + jxp) / (12.0 * dx * dy)


def velocity(psi: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    """Velocity ``u = -psi_y``, ``v = psi_x`` on all nodes (second order).

    Centered differences inside, second-order one-sided differences on the
    walls (``np.gradient`` with ``edge_order=2``).
    """
    dpsidx, dpsidy = np.gradient(psi, dx, dy, edge_order=2)
    return -dpsidy, dpsidx


def wall_vorticity(psi: np.ndarray, zeta: np.ndarray, dx: float, dy: float, bc: str,
                   bc_ns: str | None = None) -> None:
    """Fill the wall values of ``zeta`` in place from the boundary condition.

    ``psi`` vanishes on the walls. Along the western wall ``psi_yy = 0``,
    so ``zeta = psi_xx`` there.

    * free slip: ``zeta = 0`` (no tangential stress);
    * no slip: ``psi_x = 0`` on the wall; the ghost value ``psi_{-1} = psi_1``
      gives ``zeta_0 = 2 psi_1 / dx^2`` (Thom 1933).

    ``bc`` applies to the western and eastern walls, ``bc_ns`` (default
    ``bc``) to the southern and northern walls. Corner values are set to
    zero; no interior stencil uses them.
    """
    bc_ns = bc if bc_ns is None else bc_ns
    for c in (bc, bc_ns):
        if c not in ("free-slip", "no-slip"):
            raise ValueError(f"unknown boundary condition {c!r}")
    if bc == "free-slip":
        zeta[0, :] = zeta[-1, :] = 0.0
    else:
        zeta[0, :] = 2.0 * psi[1, :] / dx**2
        zeta[-1, :] = 2.0 * psi[-2, :] / dx**2
    if bc_ns == "free-slip":
        zeta[:, 0] = zeta[:, -1] = 0.0
    else:
        zeta[:, 0] = 2.0 * psi[:, 1] / dy**2
        zeta[:, -1] = 2.0 * psi[:, -2] / dy**2
    zeta[0, 0] = zeta[0, -1] = zeta[-1, 0] = zeta[-1, -1] = 0.0


class PoissonSolver:
    """Exact inverse of the five-point Laplacian with ``psi = 0`` on the walls.

    The vectors ``sin(p pi i / nx) sin(q pi j / ny)`` diagonalize the
    discrete Laplacian with eigenvalues

    ``lambda_pq = (2 cos(p pi/nx) - 2)/dx^2 + (2 cos(q pi/ny) - 2)/dy^2``,

    so ``psi = DST^-1[ DST[zeta] / lambda ]`` costs ``O(N log N)`` and is
    exact to round-off.
    """

    def __init__(self, nx: int, ny: int, dx: float, dy: float) -> None:
        p = np.arange(1, nx)
        q = np.arange(1, ny)
        lx = (2.0 * np.cos(np.pi * p / nx) - 2.0) / dx**2
        ly = (2.0 * np.cos(np.pi * q / ny) - 2.0) / dy**2
        self.eigenvalues = lx[:, None] + ly[None, :]
        self.shape = (nx - 1, ny - 1)

    def solve(self, rhs_interior: np.ndarray) -> np.ndarray:
        """Return interior ``psi`` with ``laplacian(psi) = rhs``."""
        hat = sfft.dstn(rhs_interior, type=1, workers=-1)
        return sfft.idstn(hat / self.eigenvalues, type=1, workers=-1)

    def smallest_eigenvalue_magnitude(self) -> float:
        return float(np.min(np.abs(self.eigenvalues)))
