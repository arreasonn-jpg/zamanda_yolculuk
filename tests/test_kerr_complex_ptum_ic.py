import numpy as np
from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state


def test_shapes():
    p = kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=0.5)
    for k in ("state", "PR", "PI", "PiR", "PiI", "Phi", "Pi_Phi", "psi", "Pi_psi"):
        assert k in p
    assert p["state"].phi.shape == (9, 9, 9)
    assert p["PR"].shape == (9, 9, 9)


def test_a_too_large():
    import pytest
    with pytest.raises(ValueError):
        kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=1.5)


def test_psi_activated():
    p = kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=0.5, A_fuel=3.0)
    assert np.max(p["psi"]) > 1.0


def test_a_dependence():
    p0 = kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=0.0)
    p9 = kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=0.9)
    assert not np.allclose(p0["state"].beta, p9["state"].beta)


def test_complex_scalar_perturbative():
    p = kerr_complex_ptum_state(N=9, L=20.0, M=1.0, a=0.5,
                                  Psi_R_amp=1e-4, Psi_I_amp=0.5e-4)
    assert np.max(p["PR"]) < 1e-3
    assert np.max(p["PI"]) < 1e-3
