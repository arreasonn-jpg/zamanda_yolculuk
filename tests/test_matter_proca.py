import numpy as np
import pytest
from zt005.matter_proca import proca_sources, tt_content_estimate
from zt005.bssn_state import minkowski_state

def test_vacuum_zero():
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    src = proca_sources(Z, Z, Z, Z, Z, Z, Z, s.gtilde, s.phi, h=0.1)
    assert np.allclose(src["rho"], 0.0)

def test_electric_only_positive_rho():
    """Constant E1 -> rho = 0.5 E^2."""
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    E1 = np.ones((5, 5, 5)) * 0.5
    src = proca_sources(Z, Z, Z, Z, E1, Z, Z, s.gtilde, s.phi, h=0.1)
    # rho = 0.5 * 0.25 = 0.125
    assert np.allclose(src["rho"], 0.125)

def test_anisotropic_stress():
    """Pure E along x -> T_ij is NOT isotropic."""
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    E1 = np.ones((5, 5, 5)) * 1.0
    src = proca_sources(Z, Z, Z, Z, E1, Z, Z, s.gtilde, s.phi, h=0.1)
    # S_xx = -E1^2 + 0.5*E1^2 = -0.5
    # S_yy = S_zz = +0.5
    S = src["S_ij"]
    mid = 2
    assert S[mid, mid, mid, 0, 0] < 0
    assert S[mid, mid, mid, 1, 1] > 0
    assert S[mid, mid, mid, 2, 2] > 0
    # Spatial trace = rho (4D tracelessness: -T_00 + T_ii = 0)
    # For E=(1,0,0), rho = 0.5, spatial trace should be 0.5
    tr = S[mid, mid, mid, 0, 0] + S[mid, mid, mid, 1, 1] + S[mid, mid, mid, 2, 2]
    assert abs(tr) > 0.1  # nontrivial spatial trace
    # Anisotropy: S_xx < 0 < S_yy = S_zz
    assert S[mid, mid, mid, 0, 0] < S[mid, mid, mid, 1, 1]

def test_tt_content_nonzero():
    """Vector field has nonzero TT content (unlike scalar)."""
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    E1 = np.ones((5, 5, 5)) * 1.0
    src = proca_sources(Z, Z, Z, Z, E1, Z, Z, s.gtilde, s.phi, h=0.1)
    tt = tt_content_estimate(src["S_ij"], s.gtilde, None, None, None)
    assert np.max(tt) > 0.1

def test_mass_term_positive_rho():
    """Massive Proca with nonzero A -> positive rho contribution."""
    s = minkowski_state(5)
    Z = np.zeros((5, 5, 5))
    A1 = np.ones((5, 5, 5)) * 0.3
    src = proca_sources(Z, A1, Z, Z, Z, Z, Z, s.gtilde, s.phi, h=0.1, m2=1.0)
    # rho = 0.5 * 1.0 * 0.09 = 0.045
    assert np.allclose(src["rho"], 0.045)

def test_shapes():
    s = minkowski_state(7)
    Z = np.zeros((7, 7, 7))
    src = proca_sources(Z, Z, Z, Z, Z, Z, Z, s.gtilde, s.phi, h=0.1)
    assert src["rho"].shape == (7, 7, 7)
    assert src["S_i"].shape == (7, 7, 7, 3)
    assert src["S_ij"].shape == (7, 7, 7, 3, 3)
