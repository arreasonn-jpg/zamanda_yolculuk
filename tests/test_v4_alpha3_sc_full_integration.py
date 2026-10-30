"""Test v4.0-alpha3: BSSN + SC metric + SC condensate matter integration.

This test validates the full coupling of:
- SC regular core metric (de Sitter interior + Schwarzschild exterior)
- SC condensate scalar field (Gaussian profile)
- BSSN evolution with matter sources (bssn_rhs_matter)
- RK4 time integration (rk4_step_matter)

Expected: Hamiltonian constraint growth = 1.0000 (stable)
"""
import numpy as np
import pytest
from zt005.sc_full_ic import sc_full_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full


def test_v4_alpha3_sc_full_integration_50steps():
    """50-step evolution: H growth = 1.0000."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0
    dt = 5e-4
    steps = 50
    h = L / (N - 1)

    s, phi_f, Pi_f, V_func = sc_full_state(
        N=N, L=L, M=M, rc=rc, phi_amp_target=1e-4)

    # Mask exterior region (r > 3*rc)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = R > 3.0 * rc

    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))

    for step in range(1, steps + 1):
        s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt, V_func=V_func)

    H_final = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
    nH = float(np.mean(np.abs(H_final[mask])))
    growth = nH / H0

    # Verify stability
    assert growth == pytest.approx(1.0, abs=1e-3), f"H growth = {growth:.4f}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(phi_f))
    assert np.all(np.isfinite(Pi_f))
    assert np.all(s.alpha > 0.0)
    assert np.max(np.abs(phi_f)) < 1e-3  # amplitude stays small


def test_v4_alpha3_sc_full_integration_100steps():
    """100-step evolution: H growth = 1.0000 (longer run)."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0
    dt = 5e-4
    steps = 100
    h = L / (N - 1)

    s, phi_f, Pi_f, V_func = sc_full_state(
        N=N, L=L, M=M, rc=rc, phi_amp_target=1e-4)

    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    mask = R > 3.0 * rc

    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))

    for step in range(1, steps + 1):
        s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt, V_func=V_func)

    H_final = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
    nH = float(np.mean(np.abs(H_final[mask])))
    growth = nH / H0

    # Verify stability
    assert growth == pytest.approx(1.0, abs=1e-3), f"H growth = {growth:.4f}"
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(phi_f))
    assert np.all(np.isfinite(Pi_f))
    assert np.all(s.alpha > 0.0)
    assert np.max(np.abs(phi_f)) < 1e-3


def test_v4_alpha3_sc_full_integration_shapes():
    """Verify all arrays have correct shapes."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0

    s, phi_f, Pi_f, V_func = sc_full_state(N=N, L=L, M=M, rc=rc)

    assert s.phi.shape == (N, N, N)
    assert s.gtilde.shape == (N, N, N, 3, 3)
    assert s.K.shape == (N, N, N)
    assert s.Atilde.shape == (N, N, N, 3, 3)
    assert s.Gamma.shape == (N, N, N, 3)
    assert s.alpha.shape == (N, N, N)
    assert s.beta.shape == (N, N, N, 3)
    assert phi_f.shape == (N, N, N)
    assert Pi_f.shape == (N, N, N)


def test_v4_alpha3_sc_full_integration_metric_regular():
    """Verify metric is regular (finite Ricci, bounded alpha)."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0

    s, phi_f, Pi_f, V_func = sc_full_state(N=N, L=L, M=M, rc=rc)
    h = L / (N - 1)

    from zt005.bssn_ricci_full import ricci_scalar_physical
    R = ricci_scalar_physical(s.phi, s.gtilde, h)
    assert np.all(np.isfinite(R))

    assert np.all(s.alpha > 0.0)
    assert np.all(s.alpha <= 1.0 + 1e-10)


def test_v4_alpha3_sc_full_integration_pi_zero_init():
    """Verify Pi_field starts at zero."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0

    s, phi_f, Pi_f, V_func = sc_full_state(N=N, L=L, M=M, rc=rc)
    assert np.allclose(Pi_f, 0.0)


def test_v4_alpha3_sc_full_integration_phi_amplitude():
    """Verify scalar field amplitude is at target."""
    N = 9
    L = 20.0
    M = 1.0
    rc = 2.0
    target = 1e-4

    s, phi_f, Pi_f, V_func = sc_full_state(
        N=N, L=L, M=M, rc=rc, phi_amp_target=target)

    max_phi = np.max(np.abs(phi_f))
    assert max_phi == pytest.approx(target, rel=1e-2)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])