import numpy as np
import pytest
from zt005.bssn_ic_bowen_york import bowen_york_binary

def test_shapes():
    s = bowen_york_binary(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    assert s.phi.shape == (11, 11, 11)
    assert s.Atilde.shape == (11, 11, 11, 3, 3)

def test_atilde_nonzero():
    """Bowen-York has nonzero extrinsic curvature."""
    s = bowen_york_binary(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    assert np.max(np.abs(s.Atilde)) > 1e-6

def test_atilde_traceless():
    """Atilde should be traceless."""
    s = bowen_york_binary(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    tr = np.einsum('...ij,...ij->...', s.gtilde, s.Atilde)
    interior = tr[3:-3, 3:-3, 3:-3]
    assert np.max(np.abs(interior)) < 1e-10

def test_K_zero():
    """K = 0 (maximal slicing)."""
    s = bowen_york_binary(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    assert np.allclose(s.K, 0.0)

def test_psi_matches_BL():
    """Conformal factor same as Brill-Lindquist (leading order)."""
    from zt005.bssn_ic_binary import brill_lindquist_state
    s_by = bowen_york_binary(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    s_bl = brill_lindquist_state(N=11, L=30.0, m1=0.5, m2=0.5, d=6.0)
    assert np.allclose(s_by.phi, s_bl.phi, atol=1e-10)
