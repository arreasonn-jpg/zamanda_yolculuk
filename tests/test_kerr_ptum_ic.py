import numpy as np
import pytest
from zt005.kerr_ptum_ic import kerr_ptum_state, kerr_ptum_diagnostics


def test_ic_shapes():
    s, Phi, psi, P1, P2 = kerr_ptum_state(N=9, L=20.0, M=1.0, a=0.5)
    assert s.phi.shape == (9,9,9)
    assert Phi.shape == (9,9,9)
    assert psi.shape == (9,9,9)


def test_ic_a_too_large_raises():
    with pytest.raises(ValueError):
        kerr_ptum_state(N=9, L=20.0, M=1.0, a=1.5)


def test_ic_a0_phi_positive():
    s, Phi, psi, P1, P2 = kerr_ptum_state(N=9, L=20.0, M=1.0, a=0.0)
    assert np.all(s.phi >= 0)
    assert np.all(s.alpha > 0)


def test_ic_diagnostics_structure():
    s, Phi, psi, P1, P2 = kerr_ptum_state(N=9, L=20.0, M=1.0, a=0.5,
                                            A_fuel=3.0)
    h = 20.0/8.0
    diag = kerr_ptum_diagnostics(s, Phi, P1, psi, P2, h)
    assert "rho_min" in diag
    assert "n_neg" in diag
    assert diag["n_active"] > 0  # some cells activated


def test_ic_A_fuel_determines_activation():
    """Higher A_fuel -> more activated cells."""
    s1, Phi1, psi1, _, _ = kerr_ptum_state(
        N=9, L=20.0, M=1.0, a=0.0, A_fuel=1.5, sigma=4.0)
    s2, Phi2, psi2, _, _ = kerr_ptum_state(
        N=9, L=20.0, M=1.0, a=0.0, A_fuel=3.0, sigma=4.0)
    from zt005.ptum_field import m_eff_squared
    n1 = np.sum(m_eff_squared(Phi1) < 0)
    n2 = np.sum(m_eff_squared(Phi2) < 0)
    assert n2 >= n1


def test_ic_a_dependence():
    """Kerr spin changes beta pattern (KS has radial shift even at a=0)."""
    s0, _, _, _, _ = kerr_ptum_state(N=9, L=20.0, M=1.0, a=0.0)
    s9, _, _, _, _ = kerr_ptum_state(N=9, L=20.0, M=1.0, a=0.9)
    # KS shift is nonzero even at a=0; a>0 adds azimuthal component
    assert np.max(np.abs(s0.beta)) > 0
    assert not np.allclose(s0.beta, s9.beta)
