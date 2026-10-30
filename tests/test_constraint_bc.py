"""Tests for constraint-preserving boundary conditions (v6.0-alpha3)."""
from __future__ import annotations
import numpy as np
import pytest
from src.zt005.matter.constraint_bc import (
    apply_dirichlet_bc,
    apply_neumann_bc,
    apply_sommerfeld_bc,
    apply_constraint_preserving_bc_scalar,
    apply_constraint_preserving_bc_complex,
    apply_constraint_preserving_bc_ptum,
)


def test_dirichlet_bc():
    """Dirichlet BC sets boundary to value."""
    N = 11
    f = np.ones((N, N, N)) * 5.0
    f_bc = apply_dirichlet_bc(f, 0.0)
    assert np.allclose(f_bc[0, :, :], 0.0)
    assert np.allclose(f_bc[-1, :, :], 0.0)
    assert np.allclose(f_bc[:, 0, :], 0.0)
    assert np.allclose(f_bc[:, -1, :], 0.0)
    assert np.allclose(f_bc[:, :, 0], 0.0)
    assert np.allclose(f_bc[:, :, -1], 0.0)
    # Interior unchanged
    assert np.allclose(f_bc[1:-1, 1:-1, 1:-1], 5.0)


def test_neumann_bc():
    """Neumann BC copies interior to boundary (zero normal derivative)."""
    N = 11
    x = np.linspace(0, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    f = X**2  # f = x^2, derivative = 2x
    f_bc = apply_neumann_bc(f, h=1.0)
    # At x=0 boundary, f should equal f at x=1
    assert np.allclose(f_bc[0, :, :], f[1, :, :])
    # At x=N-1 boundary, f should equal f at x=N-2
    assert np.allclose(f_bc[-1, :, :], f[-2, :, :])


def test_sommerfeld_bc():
    """Sommerfeld BC applies outgoing wave condition."""
    N = 11
    f = np.ones((N, N, N))
    h = 0.1
    dt = 0.01
    f_bc = apply_sommerfeld_bc(f, h, dt, speed=1.0)
    # For constant field, derivative is zero, so field unchanged
    assert np.allclose(f_bc, f)


def test_constraint_preserving_scalar_dirichlet():
    """Scalar field: Dirichlet for psi, Neumann for Pi."""
    N = 11
    psi = np.ones((N, N, N)) * 0.5
    Pi = np.ones((N, N, N)) * 0.3
    h = 0.1
    dt = 0.01
    
    psi_bc, Pi_bc = apply_constraint_preserving_bc_scalar(psi, Pi, h, dt, 'dirichlet')
    
    # psi = 0 at boundary
    assert np.allclose(psi_bc[0, :, :], 0.0)
    assert np.allclose(psi_bc[-1, :, :], 0.0)
    # Pi: Neumann (copy from interior)
    assert np.allclose(Pi_bc[0, :, :], Pi[1, :, :])
    assert np.allclose(Pi_bc[-1, :, :], Pi[-2, :, :])


def test_constraint_preserving_scalar_sommerfeld():
    """Scalar field: Sommerfeld for both."""
    N = 11
    psi = np.ones((N, N, N)) * 0.5
    Pi = np.ones((N, N, N)) * 0.3
    h = 0.1
    dt = 0.01
    
    psi_bc, Pi_bc = apply_constraint_preserving_bc_scalar(psi, Pi, h, dt, 'sommerfeld')
    
    # For constant fields, should be unchanged
    assert np.allclose(psi_bc, psi)
    assert np.allclose(Pi_bc, Pi)


def test_constraint_preserving_complex():
    """Complex scalar: both components get same BC."""
    N = 11
    psi_R = np.ones((N, N, N)) * 0.5
    psi_I = np.ones((N, N, N)) * 0.3
    Pi_R = np.ones((N, N, N)) * 0.2
    Pi_I = np.ones((N, N, N)) * 0.1
    h = 0.1
    dt = 0.01
    
    psi_R_bc, psi_I_bc, Pi_R_bc, Pi_I_bc = apply_constraint_preserving_bc_complex(
        psi_R, psi_I, Pi_R, Pi_I, h, dt, 'dirichlet')
    
    # Real part: Dirichlet
    assert np.allclose(psi_R_bc[0, :, :], 0.0)
    # Imag part: Dirichlet
    assert np.allclose(psi_I_bc[0, :, :], 0.0)
    # Momenta: Neumann
    assert np.allclose(Pi_R_bc[0, :, :], Pi_R[1, :, :])
    assert np.allclose(Pi_I_bc[0, :, :], Pi_I[1, :, :])


def test_constraint_preserving_ptum():
    """PTUM fields: fuel + proxy both get BC."""
    N = 11
    Phi = np.ones((N, N, N)) * 0.5
    Pi_Phi = np.ones((N, N, N)) * 0.2
    psi = np.ones((N, N, N)) * 0.3
    Pi_psi = np.ones((N, N, N)) * 0.1
    h = 0.1
    dt = 0.01
    
    Phi_bc, Pi_Phi_bc, psi_bc, Pi_psi_bc = apply_constraint_preserving_bc_ptum(
        Phi, Pi_Phi, psi, Pi_psi, h, dt, 'dirichlet')
    
    # Fuel field: Dirichlet
    assert np.allclose(Phi_bc[0, :, :], 0.0)
    # Proxy field: Dirichlet
    assert np.allclose(psi_bc[0, :, :], 0.0)
    # Momenta: Neumann
    assert np.allclose(Pi_Phi_bc[0, :, :], Pi_Phi[1, :, :])
    assert np.allclose(Pi_psi_bc[0, :, :], Pi_psi[1, :, :])


def test_invalid_bc_type():
    """Invalid BC type should raise ValueError."""
    N = 11
    psi = np.ones((N, N, N))
    Pi = np.ones((N, N, N))
    h = 0.1
    dt = 0.01
    
    with pytest.raises(ValueError):
        apply_constraint_preserving_bc_scalar(psi, Pi, h, dt, 'invalid')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])