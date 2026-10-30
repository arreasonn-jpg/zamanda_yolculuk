"""Tests for complex scalar + BSSN RK4 driver (v6.0-alpha2)."""
from __future__ import annotations
import numpy as np
import pytest
from src.zt005.bssn_state import BSSNState
from src.zt005.bssn_rk4_complex import rk4_step_complex
from src.zt005.matter.complex_scalar import complex_scalar_profile


def _make_state(N, L):
    h = L / (N - 1)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    gtilde = np.eye(3)[None, None, None, ...].repeat(N, 0).repeat(N, 1).repeat(N, 2)
    phi = np.zeros((N, N, N))
    K = np.zeros((N, N, N))
    Atilde = np.zeros((N, N, N, 3, 3))
    Gamma = np.zeros((N, N, N, 3))
    alpha = np.ones((N, N, N))
    beta = np.zeros((N, N, N, 3))
    return BSSNState(phi=phi, gtilde=gtilde, K=K, Atilde=Atilde,
                     Gamma=Gamma, alpha=alpha, beta=beta), h


def test_vacuum_stationary():
    """Vacuum (zero fields) should remain stationary."""
    N = 11
    L = 4.0
    state, h = _make_state(N, L)
    dt = 0.01

    psi_R = np.zeros((N, N, N))
    psi_I = np.zeros((N, N, N))
    Pi_R = np.zeros((N, N, N))
    Pi_I = np.zeros((N, N, N))

    new_state, new_psiR, new_psiI, new_PiR, new_PiI = rk4_step_complex(
        state, psi_R, psi_I, Pi_R, Pi_I, h, dt)

    # Metric should remain flat
    assert np.allclose(new_state.phi, 0.0, atol=1e-10)
    assert np.allclose(new_state.K, 0.0, atol=1e-10)
    assert np.allclose(new_state.Atilde, 0.0, atol=1e-10)
    # Fields should remain zero
    assert np.allclose(new_psiR, 0.0)
    assert np.allclose(new_psiI, 0.0)
    assert np.allclose(new_PiR, 0.0)
    assert np.allclose(new_PiI, 0.0)


def test_small_perturbation_finite():
    """Small Gaussian perturbation should evolve finitely."""
    N = 17
    L = 6.0
    state, h = _make_state(N, L)
    dt = 0.01

    psi_R, psi_I, Pi_R, Pi_I = complex_scalar_profile(N, L, psi0=0.05, sigma=2.0)

    new_state, new_psiR, new_psiI, new_PiR, new_PiI = rk4_step_complex(
        state, psi_R, psi_I, Pi_R, Pi_I, h, dt, m2=1.0, lam=1.0)

    # All fields should remain finite
    assert np.all(np.isfinite(new_state.phi))
    assert np.all(np.isfinite(new_state.gtilde))
    assert np.all(np.isfinite(new_state.K))
    assert np.all(np.isfinite(new_state.Atilde))
    assert np.all(np.isfinite(new_state.Gamma))
    assert np.all(np.isfinite(new_psiR))
    assert np.all(np.isfinite(new_psiI))
    assert np.all(np.isfinite(new_PiR))
    assert np.all(np.isfinite(new_PiI))


def test_metric_backreaction():
    """Non-zero field should produce metric backreaction."""
    N = 17
    L = 6.0
    state, h = _make_state(N, L)
    dt = 0.01

    psi_R, psi_I, Pi_R, Pi_I = complex_scalar_profile(N, L, psi0=0.1, sigma=2.0)

    new_state, _, _, _, _ = rk4_step_complex(
        state, psi_R, psi_I, Pi_R, Pi_I, h, dt, m2=1.0, lam=1.0)

    # K should become non-zero due to matter source
    assert not np.allclose(new_state.K, 0.0, atol=1e-6)
    # Atilde should become non-zero
    assert not np.allclose(new_state.Atilde, 0.0, atol=1e-6)


def test_u1_symmetry_preserved():
    """U(1) rotation of initial data should give same metric evolution."""
    N = 17
    L = 6.0
    state, h = _make_state(N, L)
    dt = 0.01

    psi_R, psi_I, Pi_R, Pi_I = complex_scalar_profile(N, L, psi0=0.1, sigma=2.0)
    # Add some imaginary part
    psi_I = psi_R * 0.5
    Pi_I = Pi_R * 0.5

    new_state1, _, _, _, _ = rk4_step_complex(
        state, psi_R, psi_I, Pi_R, Pi_I, h, dt, m2=1.0, lam=1.0)

    # U(1) rotation by 45 degrees
    theta = np.pi / 4
    psi_R2 = psi_R * np.cos(theta) - psi_I * np.sin(theta)
    psi_I2 = psi_R * np.sin(theta) + psi_I * np.cos(theta)
    Pi_R2 = Pi_R * np.cos(theta) - Pi_I * np.sin(theta)
    Pi_I2 = Pi_R * np.sin(theta) + Pi_I * np.cos(theta)

    new_state2, _, _, _, _ = rk4_step_complex(
        state, psi_R2, psi_I2, Pi_R2, Pi_I2, h, dt, m2=1.0, lam=1.0)

    # Metric should be identical (U(1) symmetry)
    assert np.allclose(new_state1.phi, new_state2.phi)
    assert np.allclose(new_state1.K, new_state2.K)
    assert np.allclose(new_state1.Atilde, new_state2.Atilde)


def test_shapes():
    """Output shapes should match input shapes."""
    N = 17
    L = 6.0
    state, h = _make_state(N, L)
    dt = 0.01

    psi_R, psi_I, Pi_R, Pi_I = complex_scalar_profile(N, L, psi0=0.1, sigma=2.0)

    new_state, new_psiR, new_psiI, new_PiR, new_PiI = rk4_step_complex(
        state, psi_R, psi_I, Pi_R, Pi_I, h, dt)

    assert new_state.phi.shape == (N, N, N)
    assert new_state.gtilde.shape == (N, N, N, 3, 3)
    assert new_state.K.shape == (N, N, N)
    assert new_state.Atilde.shape == (N, N, N, 3, 3)
    assert new_state.Gamma.shape == (N, N, N, 3)
    assert new_psiR.shape == (N, N, N)
    assert new_psiI.shape == (N, N, N)
    assert new_PiR.shape == (N, N, N)
    assert new_PiI.shape == (N, N, N)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
