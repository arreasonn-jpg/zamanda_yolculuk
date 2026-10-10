import numpy as np
import pytest
from zt005.sc_condensate import sc_condensate_profile, sc_condensate_potential


def test_profile_shapes():
    phi, Pi, sigma, omega, A = sc_condensate_profile(
        N=7, L=20.0, M=1.0, rho_core=1e19, m_pair_MeV=10.0)
    assert phi.shape == (7,7,7)
    assert Pi.shape == (7,7,7)
    assert sigma > 0
    assert omega > 0
    assert A > 0


def test_profile_positive_center():
    phi, Pi, sigma, omega, A = sc_condensate_profile(
        N=7, L=20.0, M=1.0, rho_core=1e19, m_pair_MeV=10.0)
    c = phi[3,3,3]
    assert c > 0
    # center value > edge value
    assert c > phi[0,0,0]


def test_profile_pi_zero():
    phi, Pi, sigma, omega, A = sc_condensate_profile(
        N=7, L=20.0, M=1.0, rho_core=1e19, m_pair_MeV=10.0)
    assert np.allclose(Pi, 0.0)


def test_potential_positive():
    V = sc_condensate_potential(m_pair_MeV=10.0, g_int=0.1)
    assert V(1.0) > 0
    assert V(0.0) == 0.0
    # monotonic in phi > 0
    assert V(2.0) > V(1.0)


def test_potential_scales_with_phi_squared():
    V = sc_condensate_potential(m_pair_MeV=10.0, g_int=0.0)
    # pure mass term: V(2x) = 4 V(x)
    assert V(2.0) == pytest.approx(4.0 * V(1.0), rel=1e-6)
