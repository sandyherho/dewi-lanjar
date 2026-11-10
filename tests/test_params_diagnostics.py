"""Parameters, scales, forcing and diagnostics."""

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dewi_lanjar import (
    DimensionalScales, Grid, Params, WindForcing, beta_at_latitude, diagnostics, solve_steady,
)


def test_params_validation():
    with pytest.raises(ValueError):
        Params(delta_s=-1)
    with pytest.raises(ValueError):
        Params(bc="sticky")
    with pytest.raises(ValueError):
        Params(bc_ns="glue")
    assert Params(bc="no-slip").bc_zonal_walls == "no-slip"
    assert Params(bc="no-slip", bc_ns="free-slip").bc_zonal_walls == "free-slip"


def test_regime_labels():
    assert "Stommel" in Params(delta_s=0.1).regime()["label"]
    assert "Munk" in Params(delta_m=0.1).regime()["label"]
    assert "f-plane" in Params(delta_s=0.1, beta_factor=0).regime()["label"]
    assert "inertial" in Params(delta_m=0.01, delta_i=0.05).regime()["label"]
    assert Params(delta_m=0.02, delta_i=0.04).regime()["Re_munk"] == pytest.approx(8.0)


def test_dimensional_scales_textbook_values():
    s = DimensionalScales(L=5e6, beta=2e-11, r=1e-6, A_h=1000.0, tau0=0.1, rho0=1025, H=1000)
    assert s.U == pytest.approx(0.1 / (1025 * 1000 * 2e-11 * 5e6))
    assert s.delta_s == pytest.approx(1e-6 / (2e-11 * 5e6))
    assert s.delta_m * s.L == pytest.approx((1000 / 2e-11) ** (1 / 3))
    assert s.sverdrups(1.0) == pytest.approx(0.1 / (1025 * 2e-11) / 1e6)
    assert "Sverdrup" in s.explain()


def test_beta_at_latitude():
    assert beta_at_latitude(0) == pytest.approx(2.289e-11, rel=1e-3)
    assert beta_at_latitude(90) == pytest.approx(0.0, abs=1e-25)


@settings(max_examples=30, deadline=None)
@given(prof=st.sampled_from(["single-gyre", "double-gyre", "asymmetric-double-gyre",
                             "triple-gyre", "cyclonic-gyre"]),
       Ly=st.floats(0.5, 2.0), amp=st.floats(0.1, 3.0), asym=st.floats(0, 1))
def test_wind_curl_is_minus_dtau_dy(prof, Ly, amp, asym):
    f = WindForcing(prof, amplitude=amp, asymmetry=asym, Ly=Ly)
    y = np.linspace(0, Ly, 4001)
    dtau = np.gradient(f.tau_x(y), y, edge_order=2)
    np.testing.assert_allclose(-dtau, f.curl_steady(0.3, y), atol=1e-4 * amp * 10 / Ly)


def test_seasonal_modulation():
    f = WindForcing(seasonal_amplitude=0.5, seasonal_period=40)
    assert f.modulation(10) == pytest.approx(1.5)
    assert f.modulation(30) == pytest.approx(0.5)


def test_normalized_entropy_limits():
    assert diagnostics.normalized_entropy(np.ones(100)) == pytest.approx(1.0)
    w = np.zeros(100)
    w[3] = 1
    assert diagnostics.normalized_entropy(w) == 0.0


def test_diagnostics_steady_balances():
    g = Grid(128, 128)
    p = Params(delta_s=0.02, delta_m=0.02, bc="no-slip")
    f = WindForcing("double-gyre")
    s = solve_steady(g, p, f)
    d = diagnostics.compute(s.psi, s.zeta, g, p, f)
    assert abs(d["energy_imbalance"]) < 1e-9          # W = dissipation
    vb = d["vorticity_budget"]
    for region in ("basin", "gyre", "west", "interior"):
        assert abs(vb[region]["tendency"]) < 1e-9 * abs(vb["gyre"]["wind"])
    # the interior is in Sverdrup balance to within interior friction ...
    assert vb["interior"]["beta"] / -vb["interior"]["wind"] > 0.6
    # ... and in the western strip beta is balanced by friction
    fr = vb["west"]["drag"] + vb["west"]["viscosity"]
    assert fr / -vb["west"]["beta"] > 0.9
    assert d["gyre_band"][1] <= 0.5 + 1e-12          # double gyre: band stops at mid-basin
    assert d["points_per_layer"] > 3


def test_theory_vs_measured_width_munk():
    g = Grid(256, 64)
    p = Params(delta_m=0.03, bc="no-slip", bc_ns="free-slip")
    s = solve_steady(g, p, WindForcing())
    d = diagnostics.compute(s.psi, s.zeta, g, p, s.forcing)
    assert d["wbc_zero_crossing"] == pytest.approx(d["theory_zero_crossing"], rel=0.06)
    assert d["wbc_e_fold"] == pytest.approx(d["theory_e_fold"], rel=0.08)


def test_describe_runs_with_and_without_scales():
    g = Grid(32, 32)
    p = Params(delta_s=0.05)
    s = solve_steady(g, p, WindForcing())
    d = diagnostics.compute(s.psi, s.zeta, g, p, s.forcing)
    txt = diagnostics.describe(d)
    assert "Energy budget" in txt
    txt2 = diagnostics.describe(d, DimensionalScales())
    assert "Sv" in txt2 and "km" in txt2
    assert math.isfinite(d["KE"])
