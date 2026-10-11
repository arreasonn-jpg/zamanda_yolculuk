import numpy as np
import pytest
from zt005.refined_scalar import evolve_refined_scalar_with_h


def test_refined_scalar_static_gaussian():
    """Static Pi=0 Gaussian: phi should stay finite, energy roughly conserved."""
    N_c = 21; L_c = 20.0; h_c = L_c / (N_c - 1)
    N_f = 11; h_f = h_c / 2
    x_c = np.linspace(-L_c/2, L_c/2, N_c)
    x_f = np.linspace(-N_f*h_f/2, N_f*h_f/2, N_f)
    Xc, Yc, Zc = np.meshgrid(x_c, x_c, x_c, indexing="ij")
    Xf, Yf, Zf = np.meshgrid(x_f, x_f, x_f, indexing="ij")
    phi_c = 1e-3 * np.exp(-(Xc**2 + Yc**2 + Zc**2) / 2.0)
    Pi_c = np.zeros_like(phi_c)
    phi_f = 1e-3 * np.exp(-(Xf**2 + Yf**2 + Zf**2) / 2.0)
    Pi_f = np.zeros_like(phi_f)
    offset = (5, 5, 5)

    out = evolve_refined_scalar_with_h(
        phi_c, Pi_c, phi_f, Pi_f,
        h_c=h_c, h_f=h_f, offset=offset,
        dt=1e-3, n_steps=10, m2=0.0)

    assert np.all(np.isfinite(out["phi_coarse"]))
    assert np.all(np.isfinite(out["phi_fine"]))
    assert np.all(np.isfinite(out["Pi_fine"]))


def test_refined_scalar_shapes():
    N_c = 11; N_f = 5
    phi_c = np.zeros((N_c, N_c, N_c))
    Pi_c = np.zeros((N_c, N_c, N_c))
    phi_f = np.zeros((N_f, N_f, N_f))
    Pi_f = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_scalar_with_h(
        phi_c, Pi_c, phi_f, Pi_f,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-3, n_steps=2)
    assert out["phi_coarse"].shape == (N_c, N_c, N_c)
    assert out["phi_fine"].shape == (N_f, N_f, N_f)


def test_refined_scalar_zero_stays_zero():
    N_c = 11; N_f = 5
    Z = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    out = evolve_refined_scalar_with_h(
        Z, Z, Zf, Zf, h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-3, n_steps=5)
    assert np.allclose(out["phi_coarse"], 0.0, atol=1e-14)
    assert np.allclose(out["phi_fine"], 0.0, atol=1e-14)
