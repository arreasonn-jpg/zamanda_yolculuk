import numpy as np
import pytest
from zt005.bssn import BSSN
from zt005.bssn_rhs_full import (
    dt_gtilde, dt_phi, dt_K, dt_Atilde, bssn_rhs_full,
)


def test_dt_gtilde_zero_when_Atilde_zero():
    assert np.allclose(dt_gtilde(1.0, np.zeros((3, 3))), 0.0)


def test_dt_phi_zero_when_K_zero():
    assert dt_phi(1.0, 0.0) == 0.0


def test_dt_K_zero_on_minkowski():
    R = np.zeros((3, 3))
    A = np.zeros((3, 3))
    assert dt_K(1.0, A, R, 0.0, 0.0) == 0.0


def test_dt_Atilde_tracefree():
    A = np.diag([0.1, -0.05, -0.05])
    R = np.diag([0.01, 0.02, -0.03])
    dAt = dt_Atilde(1.0, 0.0, A, R, 0.0)
    assert np.trace(dAt) == pytest.approx(0.0, abs=1e-12)


def test_bssn_rhs_full_keys():
    b = BSSN(0.0, np.eye(3), 0.0, np.zeros((3, 3)), np.zeros(3))
    rhs = bssn_rhs_full(b, 0.0, np.zeros(3),
                        dgtilde=np.zeros((3, 3, 3)),
                        dG=np.zeros((3, 3, 3, 3)))
    for k in ("dphi", "dgt", "dK", "dAt", "dGamma"):
        assert k in rhs


def test_bssn_rhs_full_minkowski_stationary():
    b = BSSN(0.0, np.eye(3), 0.0, np.zeros((3, 3)), np.zeros(3))
    rhs = bssn_rhs_full(b, 0.0, np.zeros(3),
                        dgtilde=np.zeros((3, 3, 3)),
                        dG=np.zeros((3, 3, 3, 3)))
    assert rhs["dphi"] == pytest.approx(0.0)
    assert np.allclose(rhs["dgt"], 0.0)
    assert rhs["dK"] == pytest.approx(0.0)
    assert np.allclose(rhs["dAt"], 0.0)
