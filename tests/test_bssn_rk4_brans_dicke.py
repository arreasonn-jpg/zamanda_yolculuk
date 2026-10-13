import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_brans_dicke import rk4_step_brans_dicke


def test_vacuum_stationary():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    s2, ph, Pph = rk4_step_brans_dicke(s, Z + 1.0, Z, h=0.1, dt=1e-5)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(ph, 1.0, atol=1e-12)


def test_pi_perturbation_finite():
    N = 7; L = 20.0; h = L/(N-1)
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z2**2
    phi_BD = 1.0 + 1e-4*np.exp(-R2/2.0)
    Pi_BD = 1e-4*np.exp(-R2/2.0)
    s2, ph, Pph = rk4_step_brans_dicke(s, phi_BD, Pi_BD, h=h, dt=1e-5)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(ph))
    assert np.all(np.isfinite(Pph))


def test_returns_three_tuple():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    out = rk4_step_brans_dicke(s, Z + 1.0, Z, h=0.1, dt=1e-5)
    assert len(out) == 3


def test_gr_limit_large_omega():
    """With omega -> infinity, BD source vanishes."""
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    s2, ph, Pph = rk4_step_brans_dicke(s, Z + 1.0, Z, h=0.1, dt=1e-5,
                                        omega=1e10)
    assert np.allclose(s2.K, 0.0, atol=1e-12)
