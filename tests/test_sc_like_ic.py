import numpy as np
import pytest
from zt005.sc_like_ic import sc_like_state, mass_function_profile
from zt005.bssn_ricci_full import ricci_scalar_physical


def test_phi_finite_at_origin():
    s = sc_like_state(N=9, L=20.0, M=1.0, rc=1.0)
    assert np.all(np.isfinite(s.phi))


def test_alpha_bounded():
    s = sc_like_state(N=9, L=20.0, M=1.0, rc=1.0)
    assert np.all(s.alpha > 0.0)
    assert np.all(s.alpha <= 1.0 + 1e-10)


def test_center_regular():
    s = sc_like_state(N=11, L=20.0, M=1.0, rc=1.0)
    # phi at center vs just outside
    c = s.phi[5,5,5]
    # center value ~ 0 in smoothed core (M_eff(0)=0)
    assert abs(c) < 0.1


def test_far_field_matches_schwarzschild():
    N = 21; L = 20.0; M = 1.0; rc = 1.0
    s = sc_like_state(N, L, M=M, rc=rc)
    # At corner r ~ sqrt(3)*10
    r_corner = np.sqrt(3.0) * 10.0
    phi_expected = np.log(1.0 + M/(2.0*r_corner))
    assert s.phi[0,0,0] == pytest.approx(phi_expected, rel=1e-4)


def test_ricci_finite():
    s = sc_like_state(N=11, L=20.0, M=1.0, rc=2.0)
    h = 20.0/(11-1)
    R = ricci_scalar_physical(s.phi, s.gtilde, h)
    assert np.all(np.isfinite(R))
    # Interior away from edges
    assert np.max(np.abs(R[3:-3,3:-3,3:-3])) < 1.0


def test_mass_function_limits():
    R, M_eff = mass_function_profile(N=51, L=20.0, M=1.0, rc=1.0)
    # Center: M_eff ~ 0
    assert M_eff[25,25,25] < 0.1
    # Far: M_eff ~ M
    assert M_eff[0,0,0] == pytest.approx(1.0, abs=0.01)
