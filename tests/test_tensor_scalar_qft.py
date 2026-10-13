import numpy as np
import pytest
from zt005.tensor_scalar_qft import (
    G_eff_factor, qft_window_prediction,
    ptum_rho_min_scaling, ford_roman_ratio_from_rho,
)


def test_G_eff_ratio_GR_limit():
    """F_phi = 1 -> G_eff = G."""
    assert G_eff_factor(1.0) == 1.0


def test_G_eff_ratio_BD():
    """Brans-Dicke: F_phi = phi."""
    assert G_eff_factor(2.0) == 0.5


def test_qft_window_invariant():
    """QFT-OK window is F(phi)-independent."""
    w1 = qft_window_prediction(1.0)
    w2 = qft_window_prediction(2.0)
    w3 = qft_window_prediction(0.5)
    assert w1["A_min"] == w2["A_min"] == w3["A_min"]
    assert w1["A_max"] == w2["A_max"] == w3["A_max"]


def test_qft_window_range():
    w = qft_window_prediction(1.0)
    assert w["A_min"] == 1.35
    assert w["A_max"] == 1.87


def test_rho_min_scaling_monotonic():
    """rho_min decreases (more negative) with A."""
    rho_16 = ptum_rho_min_scaling(1.6)
    rho_18 = ptum_rho_min_scaling(1.8)
    assert rho_18 < rho_16 < 0


def test_ford_roman_ratio_linear():
    """ratio = 0.795 * |rho_min|."""
    assert ford_roman_ratio_from_rho(-1.0) == pytest.approx(0.795)
    assert ford_roman_ratio_from_rho(-2.0) == pytest.approx(1.59)


def test_window_boundary_ratio():
    """At A_min=1.35, rho_min ~ -1.22 -> ratio slightly < 1."""
    rho = ptum_rho_min_scaling(1.35)
    r = ford_roman_ratio_from_rho(rho)
    # Window boundary: ratio near 1
    assert 0.5 < r < 1.5


def test_all_scalar_tensor_models_same_window():
    """Test multiple F(phi) values give identical window."""
    F_values = [0.3, 0.5, 1.0, 2.0, 5.0]
    windows = [qft_window_prediction(F) for F in F_values]
    A_min_values = [w["A_min"] for w in windows]
    A_max_values = [w["A_max"] for w in windows]
    assert len(set(A_min_values)) == 1
    assert len(set(A_max_values)) == 1
