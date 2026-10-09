import numpy as np
import pytest
from zt005.bssn import BSSN
from zt005.bssn_rhs import (
    conformal_christoffel, conformal_connection, ricci_tensor,
    constraint_damping_term, bssn_rhs,
)


def test_christoffel_flat_zero():
    gt = np.eye(3)
    dg = np.zeros((3, 3, 3))
    G = conformal_christoffel(gt, dg)
    assert np.allclose(G, 0.0)


def test_connection_flat_zero():
    G = np.zeros((3, 3, 3))
    g = conformal_connection(G)
    assert np.allclose(g, 0.0)


def test_ricci_flat_zero():
    G = np.zeros((3, 3, 3))
    dG = np.zeros((3, 3, 3, 3))
    R = ricci_tensor(G, dG)
    assert np.allclose(R, 0.0)


def test_damping_zero_when_constraints_zero():
    dK, dG = constraint_damping_term(0.0, np.zeros(3), eta=0.5)
    assert dK == 0.0
    assert np.allclose(dG, 0.0)


def test_damping_signs():
    dK, dG = constraint_damping_term(1.0, np.array([1.0, 0.0, 0.0]), eta=0.5)
    assert dK == pytest.approx(-0.5)
    assert dG[0] == pytest.approx(-0.5)


def test_bssn_rhs_returns_all_keys():
    b = BSSN(phi=0.0, gtilde=np.eye(3), K=0.0,
             Atilde=np.zeros((3, 3)), Gamma=np.zeros(3))
    rhs = bssn_rhs(b, rho=0.0, S=np.zeros(3),
                   dgtilde=np.zeros((3, 3, 3)),
                   dG=np.zeros((3, 3, 3, 3)))
    for k in ("dphi", "dgt", "dK", "dAt", "dGamma"):
        assert k in rhs


def test_bssn_rhs_flat_minkowski_stationary():
    b = BSSN(phi=0.0, gtilde=np.eye(3), K=0.0,
             Atilde=np.zeros((3, 3)), Gamma=np.zeros(3))
    rhs = bssn_rhs(b, rho=0.0, S=np.zeros(3),
                   dgtilde=np.zeros((3, 3, 3)),
                   dG=np.zeros((3, 3, 3, 3)))
    assert rhs["dphi"] == pytest.approx(0.0)
    assert rhs["dK"] == pytest.approx(0.0)
    assert np.allclose(rhs["dAt"], 0.0)
