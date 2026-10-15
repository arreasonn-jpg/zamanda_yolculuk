import numpy as np
import pytest
from zt005.bssn_ic_binary import (
    brill_lindquist_state, orbital_frequency,
    gw_frequency, adm_mass,
)

def test_state_shapes():
    s = brill_lindquist_state(N=11, L=20.0, m1=0.5, m2=0.5, d=4.0)
    assert s.phi.shape == (11, 11, 11)
    assert s.alpha.shape == (11, 11, 11)

def test_two_punctures_in_psi():
    """psi should have two peaks."""
    s = brill_lindquist_state(N=21, L=20.0, m1=0.5, m2=0.5, d=6.0)
    # Slice through z=0, y=0
    mid = 21 // 2
    psi = np.exp(s.phi[:, mid, mid])
    # Should have two local maxima
    peaks = np.where((psi[1:-1] > psi[:-2]) & (psi[1:-1] > psi[2:]))[0] + 1
    assert len(peaks) >= 2

def test_orbital_frequency_scaling():
    """f_orbit ~ d^{-3/2}."""
    f1 = orbital_frequency(0.5, 0.5, 4.0)
    f2 = orbital_frequency(0.5, 0.5, 16.0)
    assert f1 / f2 == pytest.approx(8.0, rel=1e-10)

def test_gw_frequency_double():
    """f_GW = 2 f_orbit."""
    f_orbit = orbital_frequency(0.5, 0.5, 4.0)
    f_gw = gw_frequency(0.5, 0.5, 4.0)
    assert f_gw == pytest.approx(2.0 * f_orbit)

def test_adm_mass_decreases_with_d():
    """Closer punctures -> larger binding energy -> smaller M_ADM."""
    m_far = adm_mass(0.5, 0.5, 20.0)
    m_near = adm_mass(0.5, 0.5, 4.0)
    assert m_near < m_far

def test_alpha_bounded():
    s = brill_lindquist_state(N=11, L=20.0, m1=0.5, m2=0.5, d=4.0)
    assert np.all(s.alpha > 0)
    assert np.all(s.alpha <= 1.0 + 1e-10)
