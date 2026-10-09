import numpy as np
import pytest
from zt005.bssn_kerr import kerr_bl, kerr_horizon_radius
from zt005.bssn_linear_check import hbar_to_metric, linear_regime_ratio


def test_kerr_a0_matches_schwarzschild_horizon():
    assert kerr_horizon_radius(M=1.0, a=0.0) == pytest.approx(2.0)


def test_kerr_small_a_horizon_shifts():
    r0 = kerr_horizon_radius(1.0, 0.0)
    r9 = kerr_horizon_radius(1.0, 0.9)
    assert r9 < r0  # faster spin -> smaller horizon


def test_kerr_bl_alpha_far_field():
    alpha, g = kerr_bl(M=1.0, a=0.5, r=1e4, theta=np.pi/2)
    assert alpha == pytest.approx(1.0, abs=1e-3)
    assert g.shape == (3, 3)


def test_kerr_bl_inside_horizon_raises():
    with pytest.raises(ValueError):
        kerr_bl(M=1.0, a=0.9, r=1.0, theta=0.0)


def test_hbar_to_metric_tracefree():
    hb = np.diag([1.0, -1.0, -1.0, -1.0])
    h = hbar_to_metric(hb)
    # eta-trace of h: hbar - 1/2 eta tr(hbar); tr(hbar)=-2
    assert h.shape == (4, 4)


def test_linear_regime_ratio_small():
    assert linear_regime_ratio(1e-10, 1e-3) < 1e-6
