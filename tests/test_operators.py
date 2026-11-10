"""Discrete operators: exactness, conservation and convergence."""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dewi_lanjar import Grid
from dewi_lanjar.operators import (
    PoissonSolver, arakawa_jacobian, ddx, ddy, laplacian, velocity, wall_vorticity,
)


def random_state(nx, ny, Lx, Ly, seed):
    g = Grid(nx, ny, Lx, Ly)
    rng = np.random.default_rng(seed)
    psi = g.pad(rng.standard_normal(g.interior_shape))
    zeta = rng.standard_normal(g.shape)
    return g, psi, zeta


@settings(max_examples=30, deadline=None)
@given(nx=st.integers(6, 40), ny=st.integers(6, 40), Lx=st.floats(0.5, 2), Ly=st.floats(0.5, 2),
       seed=st.integers(0, 10_000))
def test_arakawa_conserves_energy(nx, ny, Lx, Ly, seed):
    """sum(psi J(psi, zeta)) = 0 whenever psi = 0 on the walls (any zeta)."""
    g, psi, zeta = random_state(nx, ny, Lx, Ly, seed)
    J = arakawa_jacobian(psi, zeta, g.dx, g.dy)
    scale = np.sum(np.abs(psi[1:-1, 1:-1] * J))
    assert abs(np.sum(psi[1:-1, 1:-1] * J)) <= 1e-12 * scale


def test_arakawa_antisymmetry_and_enstrophy_with_periodic_like_interior():
    """J(a, b) = -J(b, a) and sum(zeta J) = 0 when zeta also vanishes on the walls."""
    g, psi, zeta = random_state(30, 24, 1.0, 0.8, 3)
    zeta = g.pad(zeta[1:-1, 1:-1])
    J1 = arakawa_jacobian(psi, zeta, g.dx, g.dy)
    J2 = arakawa_jacobian(zeta, psi, g.dx, g.dy)
    np.testing.assert_allclose(J1, -J2, atol=1e-10 * np.abs(J1).max())
    scale = np.sum(np.abs(zeta[1:-1, 1:-1] * J1))
    assert abs(np.sum(zeta[1:-1, 1:-1] * J1)) <= 1e-12 * scale


def test_arakawa_second_order():
    errs = []
    for n in (32, 64, 128):
        g = Grid(n, n)
        X, Y = g.X, g.Y
        a = np.sin(np.pi * X) * np.sin(np.pi * Y)
        b = np.cos(2 * X) * np.exp(Y)
        exact = (np.pi * np.cos(np.pi * X) * np.sin(np.pi * Y) * np.cos(2 * X) * np.exp(Y)
                 - np.pi * np.sin(np.pi * X) * np.cos(np.pi * Y) * (-2 * np.sin(2 * X) * np.exp(Y)))
        J = arakawa_jacobian(a, b, g.dx, g.dy)
        errs.append(np.abs(J - exact[1:-1, 1:-1]).max())
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > 1.9)


@settings(max_examples=25, deadline=None)
@given(nx=st.integers(5, 50), ny=st.integers(5, 50), Lx=st.floats(0.3, 3), Ly=st.floats(0.3, 3),
       seed=st.integers(0, 10_000))
def test_poisson_solver_is_exact_inverse(nx, ny, Lx, Ly, seed):
    g = Grid(nx, ny, Lx, Ly)
    rng = np.random.default_rng(seed)
    rhs = rng.standard_normal(g.interior_shape)
    psi = g.pad(PoissonSolver(nx, ny, g.dx, g.dy).solve(rhs))
    lap = laplacian(psi, g.dx, g.dy)
    np.testing.assert_allclose(lap, rhs, atol=1e-9 * np.abs(rhs).max())


def test_poisson_eigenvalues_negative():
    ps = PoissonSolver(20, 30, 0.05, 0.04)
    assert np.all(ps.eigenvalues < 0)


def test_centered_differences_second_order():
    for op, d in ((ddx, lambda X, Y: np.cos(X + 2 * Y)), (ddy, lambda X, Y: 2 * np.cos(X + 2 * Y))):
        errs = []
        for n in (20, 40, 80):
            g = Grid(n, n)
            f = np.sin(g.X + 2 * g.Y)
            h = g.dx if op is ddx else g.dy
            errs.append(np.abs(op(f, h) - d(g.X, g.Y)[1:-1, 1:-1]).max())
        assert np.log2(errs[0] / errs[1]) > 1.9 and np.log2(errs[1] / errs[2]) > 1.9


def test_velocity_is_divergence_free_inside():
    g, psi, _ = random_state(30, 30, 1, 1, 9)
    u, v = velocity(psi, g.dx, g.dy)
    # centered divergence of centered velocities vanishes identically
    div = (u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * g.dx) + (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * g.dy)
    assert np.abs(div[1:-1, 1:-1]).max() < 1e-9 * np.abs(u).max() / g.dx


def test_wall_vorticity_rules():
    g, psi, zeta = random_state(12, 10, 1, 1, 1)
    z = zeta.copy()
    wall_vorticity(psi, z, g.dx, g.dy, "free-slip")
    assert np.all(z[0] == 0) and np.all(z[:, -1] == 0)
    z = zeta.copy()
    wall_vorticity(psi, z, g.dx, g.dy, "no-slip")
    np.testing.assert_allclose(z[0, 1:-1], 2 * psi[1, 1:-1] / g.dx**2)
    np.testing.assert_allclose(z[1:-1, -1], 2 * psi[1:-1, -2] / g.dy**2)
    z = zeta.copy()
    wall_vorticity(psi, z, g.dx, g.dy, "no-slip", "free-slip")
    assert np.all(z[:, 0] == 0)
    np.testing.assert_allclose(z[-1, 1:-1], 2 * psi[-2, 1:-1] / g.dx**2)
    with pytest.raises(ValueError):
        wall_vorticity(psi, z, g.dx, g.dy, "sticky")
