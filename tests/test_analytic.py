"""Exact separable solution, Sverdrup interior and boundary-layer theory."""

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dewi_lanjar import Params, WindForcing, exact_linear, theory_widths
from dewi_lanjar.analytic import (
    boundary_layer_profile, characteristic_roots, measure_widths, solve_mode, sverdrup_interior,
)

mpmath = pytest.importorskip("mpmath")


def stommel_mp(x, y, ds, Lx=1.0, Ly=1.0):
    """Independent 50-digit evaluation of Stommel's (1948) closed form."""
    mpmath.mp.dps = 50
    k = mpmath.pi / Ly
    s = mpmath.sqrt(1 + 4 * ds**2 * k**2)
    m1, m2 = (-1 + s) / (2 * ds), (-1 - s) / (2 * ds)
    Xp = 1 / (ds * k)                      # a = -k for the single gyre
    E1, E2 = mpmath.exp(m1 * Lx), mpmath.exp(m2 * Lx)
    B = (1 - E1) / (E1 - E2)
    A = -1 - B
    X = Xp * (1 + A * mpmath.exp(m1 * x) + B * mpmath.exp(m2 * x))
    return X * mpmath.sin(k * y)


@pytest.mark.parametrize("ds", [0.1, 0.03, 0.01])
def test_exact_matches_50_digit_stommel(ds):
    p = Params(delta_s=ds)
    f = WindForcing("single-gyre")
    xs = np.array([0.0, 0.001, 0.02, 0.3, 0.77, 1.0])
    ys = np.array([0.1, 0.5, 0.9])
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    psi = exact_linear(X, Y, p, f)["psi"]
    ref = np.array([[float(stommel_mp(mpmath.mpf(float(x)), mpmath.mpf(float(y)), mpmath.mpf(ds)))
                     for y in ys] for x in xs])
    np.testing.assert_allclose(psi, ref, rtol=1e-11, atol=1e-11 * np.abs(ref).max())


@settings(max_examples=40, deadline=None)
@given(ds=st.floats(0.005, 0.2), dm=st.one_of(st.just(0.0), st.floats(0.005, 0.1)), b=st.floats(0.0, 1.5),
       bc=st.sampled_from(["free-slip", "no-slip"]), n=st.integers(1, 3))
def test_mode_satisfies_ode_and_walls(ds, dm, b, bc, n):
    p = Params(delta_s=ds, delta_m=dm, beta_factor=b, bc=bc)
    k = n * math.pi
    a = -k
    m = solve_mode(k, a, p)
    x = np.linspace(0, 1, 41)
    X = [m.derivative(x, d) for d in range(5)]
    resid = p.dm3 * (X[4] - 2 * k**2 * X[2] + k**4 * X[0]) - ds * (X[2] - k**2 * X[0]) - b * X[1] + a
    scale = abs(a) + p.dm3 * np.abs(X[4]).max() + ds * np.abs(X[2]).max() + b * np.abs(X[1]).max()
    assert np.abs(resid).max() <= 1e-8 * scale
    assert abs(X[0][0]) < 1e-9 * np.abs(X[0]).max() + 1e-14
    assert abs(X[0][-1]) < 1e-9 * np.abs(X[0]).max() + 1e-14
    if dm > 0:
        d = 2 if bc == "free-slip" else 1
        tol = 1e-8 * np.abs(X[d]).max() + 1e-12
        assert abs(X[d][0]) < tol and abs(X[d][-1]) < tol


def test_characteristic_roots_reduce_to_stommel_quadratic():
    p = Params(delta_s=0.05)
    r = np.sort(characteristic_roots(math.pi, p).real)
    s = math.sqrt(1 + 4 * 0.05**2 * math.pi**2)
    np.testing.assert_allclose(r, np.sort([(-1 + s) / 0.1, (-1 - s) / 0.1]))


def test_exact_needs_free_slip_zonal_walls():
    with pytest.raises(ValueError):
        exact_linear(0.5, 0.5, Params(delta_m=0.05, bc="no-slip"), WindForcing())
    exact_linear(0.5, 0.5, Params(delta_m=0.05, bc="no-slip", bc_ns="free-slip"), WindForcing())


def test_exact_velocity_and_vorticity_consistent():
    p = Params(delta_s=0.03, delta_m=0.03, bc="no-slip", bc_ns="free-slip")
    f = WindForcing("asymmetric-double-gyre")
    x = np.linspace(0, 1, 8001)
    y = np.full_like(x, 0.37)
    out = exact_linear(x, y, p, f, ("psi", "v", "zeta", "u"))
    np.testing.assert_allclose(np.gradient(out["psi"], x, edge_order=2)[5:-5], out["v"][5:-5],
                               rtol=0, atol=1e-4 * np.abs(out["v"]).max())
    # no slip on the western and eastern walls
    assert abs(out["v"][0]) < 1e-8 * np.abs(out["v"]).max()
    assert abs(out["v"][-1]) < 1e-8 * np.abs(out["v"]).max()


def test_sverdrup_interior_balance_and_east_wall():
    p = Params(delta_s=0.01)
    f = WindForcing("double-gyre")
    x = np.linspace(0, 1, 201)
    y = 0.31
    psi = sverdrup_interior(x, np.full_like(x, y), p, f)
    assert abs(psi[-1]) < 1e-14
    np.testing.assert_allclose(np.gradient(psi, x), f.curl_steady(x, y), rtol=1e-10)


def test_sverdrup_interior_custom_curl_quadrature():
    f = WindForcing("custom", curl=lambda x, y: -np.pi * np.sin(np.pi * y) * (1 + 0 * x))
    g = WindForcing("single-gyre")
    p = Params(delta_s=0.02)
    x = np.linspace(0, 1, 11)
    np.testing.assert_allclose(sverdrup_interior(x, 0.4, p, f), sverdrup_interior(x, 0.4, p, g),
                               rtol=1e-12, atol=1e-12)


def test_munk_closed_forms():
    dm = 0.04
    w = theory_widths(Params(delta_m=dm, bc="no-slip"))
    assert abs(w["zero_crossing"] - 2 * math.pi * dm / math.sqrt(3)) < 1e-6
    w = theory_widths(Params(delta_m=dm, bc="free-slip"))
    assert abs(w["zero_crossing"] - 4 * math.pi * dm / (3 * math.sqrt(3))) < 1e-6
    x = np.linspace(0, 0.5, 501)
    xi = x / dm
    G, _ = boundary_layer_profile(x, Params(delta_m=dm, bc="no-slip"))
    G_ref = 1 - np.exp(-xi / 2) * (np.cos(np.sqrt(3) * xi / 2) + np.sin(np.sqrt(3) * xi / 2) / np.sqrt(3))
    np.testing.assert_allclose(G, G_ref, atol=1e-12)


def test_stommel_closed_form():
    w = theory_widths(Params(delta_s=0.07))
    assert abs(w["e_fold"] - 0.07) < 1e-5
    assert math.isnan(w["zero_crossing"])


def test_measure_widths_on_exponential():
    x = np.linspace(0, 1, 100001)
    w = measure_widths(x, np.exp(-x / 0.05) - 0.01)
    assert abs(w["x_peak"]) < 1e-12
    assert w["zero_crossing"] == pytest.approx(0.05 * math.log(100), rel=1e-6)
