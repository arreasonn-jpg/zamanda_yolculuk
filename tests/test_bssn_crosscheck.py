import numpy as np
import pytest
from zt005.bssn_crosscheck import (
    metric_perturbation_from_phi, weak_field_phi_from_mass,
    linear_regime_consistent,
)


def test_h_from_phi_scales():
    assert metric_perturbation_from_phi(0.01) == pytest.approx(0.04)


def test_phi_weak_field_small_M_over_r():
    phi = weak_field_phi_from_mass(M=1.0, r=100.0)
    # ~ M/(2r) = 0.005
    assert phi == pytest.approx(0.005, rel=1e-2)


def test_linear_regime_consistency_at_large_r():
    out = linear_regime_consistent(M=1.0, r=1e4)
    assert out["consistent"] is True
    assert out["rel_diff"] < 0.05


def test_linear_regime_inconsistent_at_small_r():
    # very close: exact phi ~ ln(1 + M/2r) differs from 2M/r
    out = linear_regime_consistent(M=1.0, r=1.0, tol_rel=0.05)
    assert out["consistent"] is False
