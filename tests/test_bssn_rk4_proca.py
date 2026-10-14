import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_proca import rk4_step_proca

def test_vacuum_stationary():
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    out = rk4_step_proca(s, Z, Z, Z, Z, Z, Z, h=0.1, dt=1e-5)
    s2 = out[0]
    assert np.allclose(s2.phi, 0.0, atol=1e-12)

def test_returns_seven_tuple():
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    out = rk4_step_proca(s, Z, Z, Z, Z, Z, Z, h=0.1, dt=1e-5)
    assert len(out) == 7

def test_wave_propagation_finite():
    N = 7; L = 20.0; h = L/(N-1)
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z2**2
    # Plane wave in E1
    E1 = 1e-4 * np.exp(-R2/2.0)
    Z = np.zeros((N, N, N))
    out = rk4_step_proca(s, Z, Z, Z, E1, Z, Z, h=h, dt=1e-5)
    for arr in out:
        if hasattr(arr, "phi"):
            assert np.all(np.isfinite(arr.phi))
        else:
            assert np.all(np.isfinite(arr))
