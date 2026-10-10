import numpy as np
from zt005.bssn_ic import schwarzschild_isotropic_state
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def test_shapes():
    s = schwarzschild_isotropic_state(N=9, L=20.0, M=1.0)
    assert s.phi.shape == (9,9,9)
    assert s.alpha.shape == (9,9,9)


def test_phi_positive():
    s = schwarzschild_isotropic_state(N=9, L=20.0, M=1.0)
    assert np.all(s.phi > 0.0)


def test_H_small_near_center():
    """H should be small if IC is constraint-satisfying."""
    s = schwarzschild_isotropic_state(N=15, L=40.0, M=1.0)
    H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h=40.0/14.0)
    # check interior away from coarse edges
    interior = H[3:-3, 3:-3, 3:-3]
    assert np.max(np.abs(interior)) < 0.5, f"max H = {np.max(np.abs(interior))}"
