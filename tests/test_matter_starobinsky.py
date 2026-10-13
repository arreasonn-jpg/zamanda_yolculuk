import numpy as np
import pytest
from zt005.matter_starobinsky import (
    V_starobinsky, dV_starobinsky, R_from_phi, slow_roll_observables,
)


def test_V_zero_at_origin():
    assert V_starobinsky(0.0) == 0.0


def test_V_plateau_large_phi():
    """V -> 1/(8 alpha) for large phi."""
    V = V_starobinsky(20.0, alpha=1.0)
    assert V == pytest.approx(1.0/8.0, abs=1e-6)


def test_dV_positive():
    """dV/dphi > 0 for phi > 0."""
    assert dV_starobinsky(1.0) > 0
    assert dV_starobinsky(3.0) > 0


def test_R_inverse():
    """R = (e^{s phi} - 1)/(2 alpha)."""
    R = R_from_phi(0.0)
    assert R == pytest.approx(0.0, abs=1e-12)


def test_slow_roll_ns():
    """Starobinsky n_s ~ 0.965 at N=55."""
    obs = slow_roll_observables(alpha=1.0, N_efolds=55)
    assert 0.960 < obs["n_s"] < 0.970


def test_slow_roll_r_small():
    """Starobinsky r ~ 0.004 at N=55."""
    obs = slow_roll_observables(alpha=1.0, N_efolds=55)
    assert obs["r"] < 0.01
