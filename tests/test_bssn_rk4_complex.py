import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex import rk4_step_complex


def test_vacuum_stationary():
    s = minkowski_state(5)
    Z = np.zeros((5,5,5))
    s2, PR, PI, PiR, PiI = rk4_step_complex(s, Z, Z, Z, Z, h=0.1, dt=1e-5)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(PR, 0.0, atol=1e-12)


def test_complex_gaussian_finite():
    N = 7; L = 20.0; h = L/(N-1)
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    PR = 1e-4 * np.exp(-R2/0.5)
    PI = 0.5e-4 * np.exp(-R2/0.5)
    PiR = np.zeros((N,N,N))
    PiI = 0.5e-4 * np.exp(-R2/0.5)
    s2, PR2, PI2, PiR2, PiI2 = rk4_step_complex(s, PR, PI, PiR, PiI, h=h, dt=1e-5)
    for arr in (PR2, PI2, PiR2, PiI2):
        assert np.all(np.isfinite(arr))
