
"""ZT-009.1 tests."""
import numpy as np
from zt005.run_zt009_1 import (
    h00_model, sensitivity_index,
    I_BASE, R_BASE, H00_BASE,
)


def test_baseline_recovers_H00():
    h = h00_model()
    assert abs(h - H00_BASE) / H00_BASE < 1e-12


def test_current_quadratic():
    h1 = h00_model(I=I_BASE)
    h2 = h00_model(I=2*I_BASE)
    assert abs(h2 / h1 - 4.0) < 1e-12


def test_radius_quadratic():
    h1 = h00_model(R=R_BASE)
    h2 = h00_model(R=2*R_BASE)
    assert abs(h2 / h1 - 4.0) < 1e-12


def test_diameter_inverse():
    h1 = h00_model(D=1.0)
    h2 = h00_model(D=2.0)
    assert abs(h1 / h2 - 2.0) < 1e-12


def test_sensitivity_current_is_2():
    """S_I = dh/dI * I/h = 2 for quadratic."""
    S = sensitivity_index("I", rel_delta=0.001)
    assert abs(S - 2.0) < 1e-3


def test_sensitivity_radius_is_2():
    S = sensitivity_index("R", rel_delta=0.001)
    assert abs(S - 2.0) < 1e-3


def test_sensitivity_diameter_is_minus_1():
    S = sensitivity_index("D", rel_delta=0.001)
    assert abs(S - (-1.0)) < 1e-3


def test_sensitivity_frequency_is_zero():
    """h does not depend on f in our model."""
    S = sensitivity_index("f", rel_delta=0.001)
    assert abs(S) < 1e-3
