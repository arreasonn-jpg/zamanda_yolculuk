import numpy as np
import pytest
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined import evolve_refined_bssn


def test_refined_bssn_minkowski_zero_stable():
    """Zero scalar, Minkowski metric: everything stays zero."""
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Z_c = np.zeros((N_c, N_c, N_c))
    Z_f = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_bssn(
        s_c, Z_c, Z_c, s_f, Z_f, Z_f,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-4, n_steps=3)
    assert np.allclose(out["state_coarse"].phi, 0.0, atol=1e-12)
    assert np.allclose(out["state_fine"].phi, 0.0, atol=1e-12)


def test_refined_bssn_small_perturbation_finite():
    """Small scalar bump on Minkowski -> all fields stay finite."""
    N_c = 11; N_f = 5
    L_c = 20.0; L_f = 10.0
    h_c = L_c / (N_c - 1); h_f = L_f / (N_f - 1)
    x_c = np.linspace(-L_c/2, L_c/2, N_c)
    x_f = np.linspace(-L_f/2, L_f/2, N_f)
    Xc, Yc, Zc = np.meshgrid(x_c, x_c, x_c, indexing="ij")
    Xf, Yf, Zf = np.meshgrid(x_f, x_f, x_f, indexing="ij")

    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    phi_c = 1e-4 * np.exp(-(Xc**2 + Yc**2 + Zc**2) / 2.0)
    phi_f = 1e-4 * np.exp(-(Xf**2 + Yf**2 + Zf**2) / 2.0)
    Pi_c = np.zeros_like(phi_c)
    Pi_f = np.zeros_like(phi_f)

    out = evolve_refined_bssn(
        s_c, phi_c, Pi_c, s_f, phi_f, Pi_f,
        h_c=h_c, h_f=h_f, offset=(3, 3, 3),
        dt=1e-5, n_steps=3)

    for key in ("state_coarse", "state_fine"):
        assert np.all(np.isfinite(out[key].phi))
        assert np.all(np.isfinite(out[key].K))
        assert np.all(np.isfinite(out[key].gtilde))
    assert np.all(np.isfinite(out["phi_coarse"]))
    assert np.all(np.isfinite(out["phi_fine"]))


def test_refined_bssn_shapes():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Z_c = np.zeros((N_c, N_c, N_c))
    Z_f = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_bssn(
        s_c, Z_c, Z_c, s_f, Z_f, Z_f,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-4, n_steps=2)
    assert out["state_coarse"].phi.shape == (N_c, N_c, N_c)
    assert out["state_fine"].phi.shape == (N_f, N_f, N_f)
    assert out["state_coarse"].gtilde.shape == (N_c, N_c, N_c, 3, 3)
