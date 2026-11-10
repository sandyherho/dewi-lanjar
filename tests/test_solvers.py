"""Steady and transient solvers: convergence, consistency and conservation."""

import math

import numpy as np
import pytest

from dewi_lanjar import (
    Grid, Params, RunConfig, TransientSolver, WindForcing, diagnostics, exact_linear, solve_steady,
)

CASES = [
    (Params(delta_s=0.05), "single-gyre"),
    (Params(delta_m=0.05), "double-gyre"),
    (Params(delta_m=0.05, bc="no-slip", bc_ns="free-slip"), "single-gyre"),
    (Params(delta_s=0.02, delta_m=0.04, bc="no-slip", bc_ns="free-slip"), "asymmetric-double-gyre"),
    (Params(delta_s=0.1, beta_factor=0.0), "triple-gyre"),
    (Params(delta_s=0.05, Ly=0.7), "single-gyre"),
]


@pytest.mark.parametrize("p,prof", CASES)
def test_steady_second_order_against_exact(p, prof):
    f = WindForcing(prof, Ly=p.Ly)
    errs = []
    for n in (32, 64, 128):
        g = Grid(n, max(8, int(round(n * p.Ly))), 1.0, p.Ly)
        s = solve_steady(g, p, f)
        assert s.residual < 1e-8
        ex = exact_linear(g.X, g.Y, p, f)["psi"]
        errs.append(np.max(np.abs(s.psi - ex)) / np.max(np.abs(ex)))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert rates[-1] > 1.85, (errs, rates)
    assert errs[-1] < 5e-3


def test_steady_requires_friction():
    with pytest.raises(ValueError):
        solve_steady(Grid(16, 16), Params(), WindForcing())


@pytest.mark.parametrize("bc", ["free-slip", "no-slip"])
def test_linear_transient_converges_to_steady_solution(bc):
    """Same operators: the time integration relaxes to the sparse solution to round-off."""
    g = Grid(40, 40)
    p = Params(delta_s=0.08, delta_m=0.04, bc=bc)
    f = WindForcing("double-gyre")
    st = solve_steady(g, p, f)
    r = TransientSolver(g, p, f).run(RunConfig(t_end=600, n_frames=2, store_dtype="float64"))
    err = np.max(np.abs(r.psi[-1] - st.psi)) / np.max(np.abs(st.psi))
    assert err < 1e-8


def test_inviscid_unforced_energy_conservation():
    """With F = 0 and no friction only the Arakawa Jacobian and beta act: KE is conserved
    up to the O(dt^4) energy error of SSP-RK3."""
    g = Grid(48, 48)
    p = Params(delta_i=0.1, beta_factor=1.0)
    f = WindForcing(amplitude=0.0)
    rng = np.random.default_rng(0)
    from scipy.ndimage import gaussian_filter
    z0 = gaussian_filter(rng.standard_normal(g.interior_shape), 2) * 50
    s = TransientSolver(g, p, f)
    r = s.run(RunConfig(t_end=5.0, n_frames=6, cfl=0.3, diag_per_frame=1), zeta0=z0)
    ke = r.series["KE"]
    assert np.max(np.abs(ke - ke[0])) / ke[0] < 1e-4


def test_energy_budget_closes_during_spinup():
    """Numerically integrated dKE/dt from the budget terms matches the KE change."""
    g = Grid(48, 48)
    p = Params(delta_s=0.05, delta_m=0.03, delta_i=0.02)
    f = WindForcing("single-gyre")
    r = TransientSolver(g, p, f).run(RunConfig(t_end=20, n_frames=2, diag_per_frame=200, cfl=0.4))
    s = r.series
    dke = np.trapezoid(s["dKE_dt"], s["t"]) if hasattr(np, "trapezoid") else np.trapz(s["dKE_dt"], s["t"])
    assert abs(dke - (s["KE"][-1] - s["KE"][0])) < 2e-3 * s["KE"][-1]


def test_f_plane_is_symmetric_and_beta_plane_is_not():
    g = Grid(64, 64)
    f = WindForcing()
    d0 = diagnostics.compute(*_steady(g, Params(delta_s=0.05, beta_factor=0.0), f), g,
                             Params(delta_s=0.05, beta_factor=0.0), f)
    d1 = diagnostics.compute(*_steady(g, Params(delta_s=0.05), f), g, Params(delta_s=0.05), f)
    assert abs(d0["ke_centroid_x"] - 0.5) < 1e-10
    assert abs(d0["intensification_index"] - 1.0) < 1e-6
    assert d1["ke_centroid_x"] < 0.2
    assert d1["intensification_index"] > 10
    assert d1["ke_zonal_entropy"] < d0["ke_zonal_entropy"]


def _steady(g, p, f):
    s = solve_steady(g, p, f)
    return s.psi, s.zeta


def test_southern_hemisphere_mirror():
    g = Grid(48, 48)
    p = Params(delta_s=0.04)
    a = solve_steady(g, p, WindForcing())
    b = solve_steady(g, p, WindForcing(sign=-1))
    np.testing.assert_allclose(b.psi, -a.psi, atol=1e-12)


def test_float_tracks_streamline_in_steady_flow():
    from dewi_lanjar import Floats
    from scipy.ndimage import map_coordinates
    g = Grid(128, 128)
    p = Params(delta_s=0.1)
    s = solve_steady(g, p, WindForcing())
    fl = Floats(g, 20, seed=2)
    psi0 = map_coordinates(s.psi, [fl.pos[:, 0] / g.dx, fl.pos[:, 1] / g.dy], order=1)
    hist = fl.advect_steady(s.psi, t_total=0.5, n_frames=50, substeps=20)
    psi1 = map_coordinates(s.psi, [hist[-1, :, 0] / g.dx, hist[-1, :, 1] / g.dy], order=1)
    assert np.max(np.abs(psi1 - psi0)) < 0.02 * np.max(np.abs(s.psi))


def test_runconfig_validation():
    with pytest.raises(ValueError):
        RunConfig(t_end=-1)
    with pytest.raises(ValueError):
        RunConfig(t_end=10, t_record=10)
    with pytest.raises(ValueError):
        RunConfig(cfl=2)


def test_blowup_is_reported():
    g = Grid(16, 16)
    s = TransientSolver(g, Params(delta_s=0.05), WindForcing())
    with pytest.raises(FloatingPointError):
        s.run(RunConfig(t_end=50, n_frames=2, dt_max=None, cfl=1.0),
              zeta0=np.full(g.interior_shape, np.nan))


@pytest.mark.slow
def test_inertial_overshoot_increases_transport():
    g = Grid(96, 96)
    f = WindForcing()
    lin = Params(delta_s=0.03, delta_m=0.03, bc="no-slip")
    non = lin.with_(delta_i=0.04)
    cfg = RunConfig(t_end=150, n_frames=2)
    d_lin = TransientSolver(g, lin, f).run(cfg).final
    d_non = TransientSolver(g, non, f).run(cfg).final
    assert d_non["psi_max"] > d_lin["psi_max"]
    assert d_non["section_y"] > d_lin["section_y"]       # gyre center moves north
    assert math.isfinite(d_non["transport_ratio"])
