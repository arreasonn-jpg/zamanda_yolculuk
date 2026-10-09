import numpy as np
import pytest
from zt005.nonlinear_gr import (
    scalar_field_stress_energy, nonlinearity_parameter,
    orders_to_nonlinear, second_order_backreaction, NonlinearBudget,
)


def test_scalar_field_vacuum_zero():
    g = np.diag([1.0, -1.0, -1.0, -1.0])
    T = scalar_field_stress_energy(1.0, np.zeros(4), g, g.copy(), V=0.0)
    assert np.allclose(T, 0.0)


def test_nonlinearity_scales_linearly():
    assert nonlinearity_parameter(np.array([1e-3])) == pytest.approx(1e-3)
    assert nonlinearity_parameter(np.array([1e-3, -2e-3])) == pytest.approx(2e-3)


def test_orders_to_nonlinear_46():
    assert orders_to_nonlinear(1e-46) == pytest.approx(46.0, abs=0.01)
    assert orders_to_nonlinear(1e-3) == pytest.approx(3.0, abs=0.01)


def test_backreaction_ratio_h_squared():
    r = second_order_backreaction(h_amp=1e-3, k_wave=1.0)
    assert r["ratio"] == pytest.approx(1e-6, rel=1e-9)


def test_budget_summary():
    s = NonlinearBudget().summary()
    assert s["orders_to_nonlinear"] == pytest.approx(46.0, abs=0.01)
