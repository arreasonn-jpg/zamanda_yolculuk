import numpy as np
import pytest
from zt005.bssn_rhs_vec import bssn_rhs_vec
from zt005.bssn_full import bssn_rhs_beta1
from zt005.bssn_state import minkowski_state
from zt005.bssn_ic_puncture import schwarzschild_puncture_state

def test_minkowski_flat_zero():
    """Minkowski: all RHS zero."""
    s = minkowski_state(5)
    h = 0.1
    rhs = bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
    for k in rhs:
        assert np.allclose(rhs[k], 0.0, atol=1e-10), f"{k} != 0"

def test_matches_loop_version():
    """Vectorized matches pointwise loop version on Schwarzschild IC."""
    N = 7; L = 40.0; h = L/(N-1)
    s = schwarzschild_puncture_state(N, L, M=1.0, eps=1.0)
    rhs_vec = bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
    rhs_loop = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
    # Compare interior (avoid boundaries where FD differs)
    for k in ("dphi", "dgt", "dK", "dAt", "dGamma"):
        v = rhs_vec[k][2:-2, 2:-2, 2:-2]
        l = rhs_loop[k][2:-2, 2:-2, 2:-2]
        assert np.allclose(v, l, atol=1e-8, rtol=1e-6), \
            f"{k}: max diff {np.max(np.abs(v-l))}"
