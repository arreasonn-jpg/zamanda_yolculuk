import numpy as np
import pytest
from zt005.bssn import (
    ADM, BSSN, adm_to_bssn, bssn_to_adm,
    hamiltonian_constraint, momentum_constraint,
    euler_step, linear_regime_check,
)


def test_minkowski_roundtrip():
    adm = ADM(alpha=1.0, beta=np.zeros(3),
              gamma=np.eye(3), K=np.zeros((3, 3)))
    b = adm_to_bssn(adm)
    assert b.phi == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(b.gtilde, np.eye(3))
    assert b.K == pytest.approx(0.0)
    assert np.allclose(b.Atilde, 0.0)
    adm2 = bssn_to_adm(b, alpha=1.0, beta=np.zeros(3))
    assert np.allclose(adm2.gamma, np.eye(3))
    assert np.allclose(adm2.K, 0.0)


def test_conformal_det_unity():
    g = np.diag([4.0, 4.0, 4.0])
    adm = ADM(1.0, np.zeros(3), g, np.zeros((3, 3)))
    b = adm_to_bssn(adm)
    assert np.linalg.det(b.gtilde) == pytest.approx(1.0, rel=1e-12)


def test_hamiltonian_minkowski_zero():
    b = BSSN(phi=0.0, gtilde=np.eye(3), K=0.0,
             Atilde=np.zeros((3, 3)), Gamma=np.zeros(3))
    assert hamiltonian_constraint(b, rho=0.0) == pytest.approx(0.0)


def test_momentum_constraint_zero_source():
    b = BSSN(phi=0.0, gtilde=np.eye(3), K=0.0,
             Atilde=np.zeros((3, 3)), Gamma=np.zeros(3))
    assert np.allclose(momentum_constraint(b, np.zeros(3)), 0.0)


def test_euler_step_phi_advances():
    b = BSSN(0.0, np.eye(3), 0.0, np.zeros((3, 3)), np.zeros(3))
    rhs = {"dphi": 1.0, "dK": 0.5}
    b2 = euler_step(b, rhs, dt=0.1)
    assert b2.phi == pytest.approx(0.1)
    assert b2.K == pytest.approx(0.05)


def test_linear_regime_threshold():
    assert linear_regime_check(1e-4) is True
    assert linear_regime_check(1.0) is False
    assert linear_regime_check(1e-46) is True
