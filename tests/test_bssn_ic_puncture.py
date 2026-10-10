import numpy as np
import pytest
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def test_smooth_no_nan():
    s = schwarzschild_puncture_state(N=11, L=20.0, M=1.0, eps=1.0)
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.alpha))


def test_phi_far_field():
    s = schwarzschild_puncture_state(N=21, L=20.0, M=1.0, eps=1.0)
    r = np.sqrt(3) * 10.0
    expected = np.log(1 + 1/(2*r))
    assert abs(s.phi[0,0,0] - expected) < 1e-3


@pytest.mark.slow
def test_H_converges_in_exterior():
    Hs = []
    for N in [21, 41, 81]:
        L = 20.0
        h = L / (N - 1)
        eps = 1.0
        s = schwarzschild_puncture_state(N, L, M=1.0, eps=eps)
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        x = np.linspace(-L/2, L/2, N)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        R = np.sqrt(X**2 + Y**2 + Z**2)
        mask = (R > 5*eps) & (R < L/3)
        Hs.append(float(np.max(np.abs(H[mask]))))
    r1 = Hs[1] / Hs[0]
    r2 = Hs[2] / Hs[1]
    print("max|H|_ext:", Hs)
    print("ratios:", r1, r2, "(expect ~4 for 2nd order)")
    assert Hs[-1] < Hs[0], "H does not decrease: " + str(Hs)
