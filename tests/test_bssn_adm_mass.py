import numpy as np
import pytest
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_adm_mass import adm_mass_approx


def test_minkowski_zero_mass():
    from zt005.bssn_state import minkowski_state
    s = minkowski_state(15)
    M = adm_mass_approx(s.phi, h=1.0)
    assert abs(M) < 1e-10


@pytest.mark.slow
def test_puncture_mass_order_of_magnitude():
    """For puncture IC, ADM mass should be ~ M within a factor of 2."""
    M_target = 1.0
    N = 41; L = 40.0; h = L/(N-1)
    s = schwarzschild_puncture_state(N, L, M=M_target, eps=1.0)
    M_adm = adm_mass_approx(s.phi, h)
    assert 0.5 * M_target < abs(M_adm) < 2.0 * M_target, f"M_adm = {M_adm}"
