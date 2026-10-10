import numpy as np
from zt005.sc_full_ic import sc_full_state
from zt005.bssn_ricci_full import ricci_scalar_physical


def test_shapes():
    s, phi_f, Pi_f, V = sc_full_state(N=9, L=20.0, M=1.0, rc=2.0)
    assert s.phi.shape == (9,9,9)
    assert phi_f.shape == (9,9,9)
    assert Pi_f.shape == (9,9,9)


def test_phi_amplitude_small():
    s, phi_f, Pi_f, V = sc_full_state(
        N=9, L=20.0, M=1.0, rc=2.0, phi_amp_target=1e-4)
    assert np.max(np.abs(phi_f)) < 1e-3


def test_pi_zero():
    s, phi_f, Pi_f, V = sc_full_state(N=9, L=20.0, M=1.0, rc=2.0)
    assert np.allclose(Pi_f, 0.0)


def test_metric_regular():
    s, phi_f, Pi_f, V = sc_full_state(N=9, L=20.0, M=1.0, rc=2.0)
    h = 20.0/(9-1)
    R = ricci_scalar_physical(s.phi, s.gtilde, h)
    assert np.all(np.isfinite(R))


def test_metric_alpha_bounded():
    s, phi_f, Pi_f, V = sc_full_state(N=9, L=20.0, M=1.0, rc=2.0)
    assert np.all(s.alpha > 0.0)
    assert np.all(s.alpha <= 1.0 + 1e-10)
