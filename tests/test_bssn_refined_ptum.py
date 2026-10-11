import numpy as np
import pytest
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_ptum import evolve_refined_ptum


def test_refined_ptum_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_ptum(
        s_c, Zc, Zc, Zc, Zc, s_f, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-4, n_steps=2)
    for key in ("state_c", "state_f"):
        assert np.allclose(out[key].phi, 0.0, atol=1e-12)


def test_refined_ptum_shapes():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_ptum(
        s_c, Zc, Zc, Zc, Zc, s_f, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-4, n_steps=1)
    assert out["state_c"].gtilde.shape == (N_c, N_c, N_c, 3, 3)
    assert out["state_f"].gtilde.shape == (N_f, N_f, N_f, 3, 3)
    assert out["Phi_c"].shape == (N_c, N_c, N_c)
    assert out["Phi_f"].shape == (N_f, N_f, N_f)


def test_refined_ptum_small_perturbation_finite():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    xc = np.linspace(-10, 10, N_c)
    xf = np.linspace(-2.5, 2.5, N_f)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    Rc = Xc**2 + Yc**2 + Zc**2
    Rf = Xf**2 + Yf**2 + Zf**2
    Phi_c = 0.5 * np.exp(-Rc/4.0)
    Phi_f = 0.5 * np.exp(-Rf/1.0)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))

    out = evolve_refined_ptum(
        s_c, Phi_c, Zc, Zc, Zc, s_f, Phi_f, Zf, Zf, Zf,
        h_c=2.0, h_f=1.25, offset=(3, 3, 3),
        dt=1e-4, n_steps=3)
    for k in ("Phi_c", "Phi_f", "psi_c", "psi_f"):
        assert np.all(np.isfinite(out[k])), f"{k} not finite"
