"""Tests for complex scalar matter module (v6.0-alpha1)."""
from __future__ import annotations
import numpy as np
import pytest
from src.zt005.matter.complex_scalar import (
    complex_scalar_sources,
    complex_scalar_field_rhs,
    complex_scalar_profile,
)


def test_vacuum_zero():
    """Zero fields -> zero sources."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.zeros((N, N, N))
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h)
    assert np.allclose(src["rho"], 0.0)
    assert np.allclose(src["S_i"], 0.0)
    assert np.allclose(src["S_ij"], 0.0)


def test_real_field_only():
    """Real field only (psi_I = 0) should match real scalar."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.ones((N, N, N)) * 0.1
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h, m2=1.0, lam=1.0)
    # V = m^2 psi^2 + lam psi^4 = 0.01 + 0.0001 = 0.0101
    # rho = V = 0.0101
    expected_rho = 1.0 * 0.1**2 + 1.0 * 0.1**4
    assert np.allclose(src["rho"], expected_rho)


def test_imaginary_field_only():
    """Imaginary field only should give same energy as real."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.zeros((N, N, N))
    psi_I = np.ones((N, N, N)) * 0.1
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h, m2=1.0, lam=1.0)
    expected_rho = 1.0 * 0.1**2 + 1.0 * 0.1**4
    assert np.allclose(src["rho"], expected_rho)


def test_momentum_contribution():
    """Non-zero Pi should contribute to rho."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.zeros((N, N, N))
    psi_I = np.zeros((N, N, N))
    Pi_R = np.ones((N, N, N)) * 0.5
    Pi_I = np.zeros((N, N, N))
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h)
    # rho = 0.5 * Pi^2 = 0.5 * 0.25 = 0.125
    assert np.allclose(src["rho"], 0.125)


def test_gradient_contribution():
    """Non-zero gradient should contribute to rho."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    psi_R = X  # linear gradient in x: f(x) = x
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h)
    # For f(x) = x, central difference gives exactly 1.0 at interior points
    # grad^2 = 1.0, rho = 0.5 * 1.0 = 0.5
    assert np.all(src["rho"] > 0)
    c = N // 2
    assert np.isclose(src["rho"][c, c, c], 0.5, rtol=1e-2)


def test_field_rhs_zero():
    """Zero fields -> zero RHS."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.zeros((N, N, N))
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    alpha = np.ones((N, N, N))
    K = np.zeros((N, N, N))

    dpsiR, dpsiI, dPiR, dPiI = complex_scalar_field_rhs(
        psi_R, psi_I, Pi_R, Pi_I, alpha, K, h)
    assert np.allclose(dpsiR, 0.0)
    assert np.allclose(dpsiI, 0.0)
    assert np.allclose(dPiR, 0.0)
    assert np.allclose(dPiI, 0.0)


def test_field_rhs_potential():
    """Potential term should drive Pi evolution."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.ones((N, N, N)) * 0.1
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))
    alpha = np.ones((N, N, N))
    K = np.zeros((N, N, N))

    dpsiR, dpsiI, dPiR, dPiI = complex_scalar_field_rhs(
        psi_R, psi_I, Pi_R, Pi_I, alpha, K, h, m2=1.0, lam=1.0)
    # dPiR/dt = -alpha * dV/dpsi_R = - (2*m^2*psi_R + 4*lam*psi^2*psi_R)
    # = - (2*1*0.1 + 4*1*0.01*0.1) = -0.204
    expected = -0.204
    assert np.allclose(dPiR, expected, atol=1e-3)


def test_profile_shapes():
    """Profile function returns correct shapes."""
    N = 17
    L = 6.0
    psi_R, psi_I, Pi_R, Pi_I = complex_scalar_profile(N, L, psi0=0.1, sigma=2.0)
    assert psi_R.shape == (N, N, N)
    assert psi_I.shape == (N, N, N)
    assert Pi_R.shape == (N, N, N)
    assert Pi_I.shape == (N, N, N)
    # Center should be psi0
    c = N // 2
    assert np.isclose(psi_R[c, c, c], 0.1, atol=1e-2)
    # Imaginary part zero
    assert np.allclose(psi_I, 0.0)
    # Momenta zero
    assert np.allclose(Pi_R, 0.0)
    assert np.allclose(Pi_I, 0.0)


def test_u1_symmetry():
    """U(1) rotation should leave rho invariant."""
    N = 11
    L = 4.0
    h = L / (N - 1)
    psi_R = np.ones((N, N, N)) * 0.1
    psi_I = np.ones((N, N, N)) * 0.2
    Pi_R = np.ones((N, N, N)) * 0.05
    Pi_I = np.ones((N, N, N)) * 0.1
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi_conf = np.zeros((N, N, N))

    src1 = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h)

    # U(1) rotation by 45 degrees
    theta = np.pi / 4
    psi_R2 = psi_R * np.cos(theta) - psi_I * np.sin(theta)
    psi_I2 = psi_R * np.sin(theta) + psi_I * np.cos(theta)
    Pi_R2 = Pi_R * np.cos(theta) - Pi_I * np.sin(theta)
    Pi_I2 = Pi_R * np.sin(theta) + Pi_I * np.cos(theta)

    src2 = complex_scalar_sources(psi_R2, psi_I2, Pi_R2, Pi_I2, gtilde, phi_conf, h)

    assert np.allclose(src1["rho"], src2["rho"])
    assert np.allclose(src1["S"], src2["S"])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])