import numpy as np
import pytest
from zt005.bssn_z4c import damping_dK, damping_dGamma, damping_energy


def test_damping_dK_zero_when_H_zero():
    assert damping_dK(1.0, 0.0) == 0.0


def test_damping_dK_sign_opposes_H():
    assert damping_dK(1.0, 1.0, kappa1=0.5) == pytest.approx(-0.5)
    assert damping_dK(1.0, -1.0, kappa1=0.5) == pytest.approx(0.5)


def test_damping_dGamma_vector():
    M = np.array([1.0, -1.0, 0.0])
    out = damping_dGamma(M, kappa2=0.1)
    assert out[0] == pytest.approx(-0.1)
    assert out[1] == pytest.approx(0.1)
    assert out[2] == 0.0


def test_damping_energy_scalar():
    assert damping_energy(0.5, np.array([0.0, 0.0, 0.0])) == pytest.approx(0.5)


def test_damping_energy_with_vector():
    e = damping_energy(0.5, np.array([0.3, 0.4, 0.0]))
    assert e == pytest.approx(0.5 + 0.5)
