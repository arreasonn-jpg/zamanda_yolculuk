import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_sc_ptum import rk4_step_sc_ptum
from zt005.matter.sc_condensate import sc_condensate_potential


def test_vacuum_stationary():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}
    out = rk4_step_sc_ptum(s, Z, Z, Z, Z, Z, Z, h=0.1, dt=1e-5,
                            V_sc=V_sc, ptum_params=params)
    s2, psc, Psc, Ph, PPh, ps, Pps = out
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(psc, 0.0, atol=1e-12)
    assert np.allclose(Ph, 0.0, atol=1e-12)


def test_returns_seven_tuple():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}
    out = rk4_step_sc_ptum(s, Z, Z, Z, Z, Z, Z, h=0.1, dt=1e-5,
                            V_sc=V_sc, ptum_params=params)
    assert len(out) == 7


def test_small_perturbation_finite():
    N = 7; L = 20.0; h = L/(N-1)
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z2**2
    phi_sc = 1e-4 * np.exp(-R2 / 2.0)
    Phi = 0.5 * np.exp(-R2 / 8.0)
    Z = np.zeros((N,N,N))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}
    out = rk4_step_sc_ptum(s, phi_sc, Z, Phi, Z, Z, Z, h=h, dt=1e-5,
                            V_sc=V_sc, ptum_params=params)
    for arr in out:
        if hasattr(arr, "phi"):
            assert np.all(np.isfinite(arr.phi))
        else:
            assert np.all(np.isfinite(arr))
