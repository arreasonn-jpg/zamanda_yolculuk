import numpy as np
from zt005.matter_ptum import ptum_sources
from zt005.bssn_state import minkowski_state
from zt005.ptum_field import fuel_profile


def test_vacuum_zero():
    s = minkowski_state(5)
    Phi = np.zeros((5,5,5))
    psi = np.zeros((5,5,5))
    Pi_Phi = np.zeros((5,5,5))
    Pi_psi = np.zeros((5,5,5))
    src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0)
    assert np.allclose(src["S_i"], 0.0)


def test_fuel_phi_only_has_positive_rho():
    """Fuel Phi = 0.5, V_fuel > 0 -> rho > 0."""
    N = 5
    s = minkowski_state(N)
    Phi = np.ones((N,N,N)) * 0.5
    psi = np.zeros((N,N,N))
    Pi_Phi = np.zeros((N,N,N))
    Pi_psi = np.zeros((N,N,N))
    src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi,
                       h=0.1, V0=1.0, lam=0.05)
    assert np.all(src["rho"] > 0)


def test_proxy_activated_negative_contribution():
    """At Phi > Phi_c, V_proxy has negative minimum at psi_0."""
    from zt005.ptum_field import psi_equilibrium
    N = 5
    s = minkowski_state(N)
    Phi = np.ones((N,N,N)) * 2.0
    psi_0 = psi_equilibrium(2.0, m2=1.0, alpha=1.0, Eth=1.0, lam_psi=1.0)
    psi = np.ones((N,N,N)) * psi_0
    Pi_Phi = np.zeros((N,N,N))
    Pi_psi = np.zeros((N,N,N))
    src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi,
                       h=0.1, V0=0.0, lam=0.05)  # turn off fuel V
    # V_proxy at minimum is negative: V = -m_eff^4/(4 lam)
    assert np.all(src["rho"] < 0), f"rho = {src['rho'].min()}"


def test_shapes():
    s = minkowski_state(7)
    Phi = fuel_profile(7, 20.0, Phi0=0.5, sigma=2.0)
    psi = np.zeros((7,7,7))
    Pi = np.zeros((7,7,7))
    src = ptum_sources(Phi, Pi, psi, Pi, s.gtilde, s.phi, h=0.1)
    assert src["rho"].shape == (7,7,7)
    assert src["S_i"].shape == (7,7,7,3)
    assert src["S_ij"].shape == (7,7,7,3,3)
