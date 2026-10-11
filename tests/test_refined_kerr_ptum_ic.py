import numpy as np
from zt005.refined_kerr_ptum_ic import refined_kerr_ptum_ic


def test_shapes():
    s_c, Phi_c, psi_c, s_f, Phi_f, psi_f = refined_kerr_ptum_ic(
        N_c=9, L_c=20.0, N_f=9, L_f=10.0, M=1.0, a=0.5)
    assert s_c.phi.shape == (9, 9, 9)
    assert s_f.phi.shape == (9, 9, 9)
    assert Phi_c.shape == (9, 9, 9)


def test_center_values_match():
    """At r=0 both grids see same A_fuel value."""
    s_c, Phi_c, psi_c, s_f, Phi_f, psi_f = refined_kerr_ptum_ic(
        N_c=9, L_c=20.0, N_f=9, L_f=10.0, M=1.0, a=0.0, A_fuel=3.0)
    # Center index (4,4,4)
    assert Phi_c[4, 4, 4] == Phi_f[4, 4, 4]


def test_a_nonzero_beta():
    s_c, _, _, s_f, _, _ = refined_kerr_ptum_ic(
        N_c=9, L_c=20.0, N_f=9, L_f=10.0, M=1.0, a=0.9)
    assert np.max(np.abs(s_c.beta)) > 0
    assert np.max(np.abs(s_f.beta)) > 0


def test_psi_equilibrium_positive_center():
    _, _, psi_c, _, _, psi_f = refined_kerr_ptum_ic(
        N_c=9, L_c=20.0, N_f=9, L_f=10.0, M=1.0, a=0.5, A_fuel=3.0)
    assert psi_c[4, 4, 4] > 0
    assert psi_f[4, 4, 4] > 0
