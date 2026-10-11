import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.ptum_field import fuel_profile


def test_vacuum_stationary():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    s2, Phi, PiPhi, psi, Pipsi = rk4_step_ptum(s, Z, Z, Z, Z, h=0.1, dt=1e-5)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(Phi, 0.0, atol=1e-12)
    assert np.allclose(psi, 0.0, atol=1e-12)


def test_fuel_profile_evolves_finite():
    N = 7; L = 20.0; h = L/(N-1)
    s = minkowski_state(N)
    Phi = fuel_profile(N, L, Phi0=0.3, sigma=2.0)
    Z = np.zeros((N,N,N))
    s2, Phi2, PiPhi2, psi2, Pipsi2 = rk4_step_ptum(
        s, Phi, Z, Z, Z, h=h, dt=1e-4)
    assert np.all(np.isfinite(Phi2))
    assert np.all(np.isfinite(psi2))
    assert np.all(np.isfinite(s2.phi))
