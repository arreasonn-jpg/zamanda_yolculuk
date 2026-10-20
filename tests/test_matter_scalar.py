import numpy as np
import pytest
from zt005.matter.matter_scalar import scalar_sources
from zt005.bssn_state import minkowski_state


def test_vacuum_zero_phi_zero_sources():
    s = minkowski_state(5)
    phi_f = np.zeros((5,5,5))
    Pi_f = np.zeros((5,5,5))
    src = scalar_sources(phi_f, Pi_f, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0)
    assert np.allclose(src["S_i"], 0.0)
    assert np.allclose(src["S_ij"], 0.0)


def test_homogeneous_pi_only():
    """Homogeneous field: Pi uniform -> rho = 1/2 Pi^2, S_ij = gamma_ij * L."""
    s = minkowski_state(5)
    phi_f = np.zeros((5,5,5))
    Pi_f = np.ones((5,5,5)) * 0.5
    src = scalar_sources(phi_f, Pi_f, s.gtilde, s.phi, h=0.1)
    # rho = 1/2 * 0.25 = 0.125
    assert np.allclose(src["rho"], 0.125)
    # S_i = 0
    assert np.allclose(src["S_i"], 0.0)
    # S = gamma^ij S_ij = 3 * (1/2 * 0.25) = 3/8
    assert np.allclose(src["S"], 0.375)


def test_homogeneous_phi_gradient():
    """Only grad phi nonzero: S_ii has anisotropic structure."""
    N = 5
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    phi_f = X.copy()
    Pi_f = np.zeros((N,N,N))
    src = scalar_sources(phi_f, Pi_f, s.gtilde, s.phi, h=2.0/(N-1))
    # |grad phi|^2 = 1
    # rho = 1/2
    assert np.allclose(src["rho"][2:-2,2:-2,2:-2], 0.5, atol=1e-10)


def test_phi_only_no_pi_isotropic_stress():
    """phi = const -> rho = V(phi), S_ij = -gamma_ij V."""
    N = 5
    s = minkowski_state(N)
    phi_f = np.ones((N,N,N)) * 0.3
    Pi_f = np.zeros((N,N,N))
    V_func = lambda p: 0.5 * p**2
    src = scalar_sources(phi_f, Pi_f, s.gtilde, s.phi, h=0.1, V_func=V_func)
    # rho = V(0.3) = 0.045
    assert np.allclose(src["rho"], 0.045)
    # S_ij = -gamma_ij V -> S = -3 V
    assert np.allclose(src["S"], -3.0 * 0.045)


def test_shapes():
    s = minkowski_state(7)
    phi_f = np.zeros((7,7,7))
    Pi_f = np.zeros((7,7,7))
    src = scalar_sources(phi_f, Pi_f, s.gtilde, s.phi, h=0.1)
    assert src["rho"].shape == (7,7,7)
    assert src["S_i"].shape == (7,7,7,3)
    assert src["S_ij"].shape == (7,7,7,3,3)
    assert src["S"].shape == (7,7,7)
