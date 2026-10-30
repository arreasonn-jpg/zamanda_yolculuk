"""Constraint-preserving boundary conditions for matter sources (v6.0-alpha3).

Implements boundary conditions that preserve the Hamiltonian and momentum
constraints at the outer boundary for matter-coupled BSSN.

For scalar fields, the constraint-preserving BCs are:
  - Dirichlet for the field (psi = 0 at boundary)
  - Neumann for the momentum (d_n Pi = 0 at boundary)
  - Sommerfeld/outgoing for radiative fields

For PTUM fields (Phi, psi):
  - Phi: Dirichlet (Phi = 0 at boundary) or Sommerfeld
  - psi: Dirichlet (psi = 0 at boundary) or Sommerfeld
  - Momenta: Neumann (d_n Pi = 0)
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central


def _normal_derivative(f, h, axis, side='outer'):
    """Normal derivative at boundary along given axis.
    side='outer': derivative at N-2 pointing outward
    side='inner': derivative at 1 pointing inward
    """
    N = f.shape[axis]
    if side == 'outer':
        # One-sided second-order at boundary
        # f'(N-1) ≈ (-3f(N-1) + 4f(N-2) - f(N-3)) / (2h)
        sl = [slice(None)] * f.ndim
        sl[axis] = -1
        sl_m1 = [slice(None)] * f.ndim
        sl_m1[axis] = -2
        sl_m2 = [slice(None)] * f.ndim
        sl_m2[axis] = -3
        return (-3.0*f[tuple(sl)] + 4.0*f[tuple(sl_m1)] - f[tuple(sl_m2)]) / (2.0*h)
    else:
        # f'(0) ≈ (3f(0) - 4f(1) + f(2)) / (2h)
        sl = [slice(None)] * f.ndim
        sl[axis] = 0
        sl_p1 = [slice(None)] * f.ndim
        sl_p1[axis] = 1
        sl_p2 = [slice(None)] * f.ndim
        sl_p2[axis] = 2
        return (3.0*f[tuple(sl)] - 4.0*f[tuple(sl_p1)] + f[tuple(sl_p2)]) / (2.0*h)


def apply_dirichlet_bc(field, value=0.0):
    """Apply Dirichlet BC: field = value at all boundaries."""
    f = field.copy()
    f[0, :, :] = value
    f[-1, :, :] = value
    f[:, 0, :] = value
    f[:, -1, :] = value
    f[:, :, 0] = value
    f[:, :, -1] = value
    return f


def apply_neumann_bc(field, h):
    """Apply Neumann BC: d_n field = 0 at all boundaries.
    Uses ghost cell method: f_ghost = f_interior
    """
    f = field.copy()
    # x boundaries
    f[0, :, :] = f[1, :, :]
    f[-1, :, :] = f[-2, :, :]
    # y boundaries
    f[:, 0, :] = f[:, 1, :]
    f[:, -1, :] = f[:, -2, :]
    # z boundaries
    f[:, :, 0] = f[:, :, 1]
    f[:, :, -1] = f[:, :, -2]
    return f


def apply_sommerfeld_bc(field, h, dt, speed=1.0):
    """Apply Sommerfeld (outgoing wave) BC:
    d_t f + speed * d_n f = 0 at boundary
    Discretized: f_new = f_old - speed * dt * d_n f
    """
    f = field.copy()
    # x boundaries
    dnx_outer = _normal_derivative(f, h, 0, 'outer')
    dnx_inner = _normal_derivative(f, h, 0, 'inner')
    f[-1, :, :] = f[-1, :, :] - speed * dt * dnx_outer
    f[0, :, :] = f[0, :, :] + speed * dt * dnx_inner
    # y boundaries
    dny_outer = _normal_derivative(f, h, 1, 'outer')
    dny_inner = _normal_derivative(f, h, 1, 'inner')
    f[:, -1, :] = f[:, -1, :] - speed * dt * dny_outer
    f[:, 0, :] = f[:, 0, :] + speed * dt * dny_inner
    # z boundaries
    dnz_outer = _normal_derivative(f, h, 2, 'outer')
    dnz_inner = _normal_derivative(f, h, 2, 'inner')
    f[:, :, -1] = f[:, :, -1] - speed * dt * dnz_outer
    f[:, :, 0] = f[:, :, 0] + speed * dt * dnz_inner
    return f


def apply_constraint_preserving_bc_scalar(psi, Pi, h, dt, bc_type='dirichlet'):
    """Constraint-preserving BC for single real scalar field.
    
    bc_type:
      'dirichlet': psi=0, d_n Pi=0
      'sommerfeld': outgoing for both
      'neumann': d_n psi=0, d_n Pi=0
    """
    if bc_type == 'dirichlet':
        psi_bc = apply_dirichlet_bc(psi, 0.0)
        Pi_bc = apply_neumann_bc(Pi, h)
    elif bc_type == 'sommerfeld':
        psi_bc = apply_sommerfeld_bc(psi, h, dt)
        Pi_bc = apply_sommerfeld_bc(Pi, h, dt)
    elif bc_type == 'neumann':
        psi_bc = apply_neumann_bc(psi, h)
        Pi_bc = apply_neumann_bc(Pi, h)
    else:
        raise ValueError(f"Unknown bc_type: {bc_type}")
    return psi_bc, Pi_bc


def apply_constraint_preserving_bc_complex(psi_R, psi_I, Pi_R, Pi_I, h, dt, bc_type='dirichlet'):
    """Constraint-preserving BC for complex scalar field."""
    psi_R_bc, Pi_R_bc = apply_constraint_preserving_bc_scalar(psi_R, Pi_R, h, dt, bc_type)
    psi_I_bc, Pi_I_bc = apply_constraint_preserving_bc_scalar(psi_I, Pi_I, h, dt, bc_type)
    return psi_R_bc, psi_I_bc, Pi_R_bc, Pi_I_bc


def apply_constraint_preserving_bc_ptum(Phi, Pi_Phi, psi, Pi_psi, h, dt, bc_type='dirichlet'):
    """Constraint-preserving BC for PTUM fields (fuel + proxy)."""
    Phi_bc, Pi_Phi_bc = apply_constraint_preserving_bc_scalar(Phi, Pi_Phi, h, dt, bc_type)
    psi_bc, Pi_psi_bc = apply_constraint_preserving_bc_scalar(psi, Pi_psi, h, dt, bc_type)
    return Phi_bc, Pi_Phi_bc, psi_bc, Pi_psi_bc


def hamiltonian_constraint_residual(phi, gtilde, K, Atilde, alpha, h, rho):
    """Compute Hamiltonian constraint residual: H = R + K^2 - K_ij K^ij - 16 pi rho.
    
    Returns the constraint violation at each point.
    """
    from .bssn_rhs import ricci_scalar_from_grid
    R = ricci_scalar_from_grid(phi, gtilde, h)
    K2 = K**2
    A2 = np.einsum('...ij,...ij->...', Atilde, Atilde) * np.exp(-4.0*phi)
    H = R + K2 - A2 - 16.0 * np.pi * rho
    return H


def momentum_constraint_residual(phi, gtilde, K, Atilde, Gamma, alpha, h, S_i):
    """Compute momentum constraint residual: M_i = D_j (A^ij - 1/3 gamma^ij K) - 8 pi S_i."""
    from .bssn_rhs import conformal_christoffel
    from .bssn_grid import dgtilde_from_grid
    
    Nx, Ny, Nz = phi.shape
    M = np.zeros((Nx, Ny, Nz, 3))
    
    e4phi = np.exp(4.0 * phi)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)
    
    # D_j A^ij = d_j A^ij + Gamma^i_jk A^kj + Gamma^j_jk A^ik - Gamma^k_jk A^ij
    # Simplified: use conformal connection functions
    # M_i = d_j Atilde^ij + 6 Atilde^ij d_j phi - 2/3 d^i K - 8 pi S_i
    
    dAtilde = np.zeros((Nx, Ny, Nz, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            comp = Atilde[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dAtilde[..., i, j, 0] = dx
            dAtilde[..., i, j, 1] = dy
            dAtilde[..., i, j, 2] = dz
    
    dphi = np.stack(d_central(phi, h), axis=-1)
    dK = np.stack(d_central(K, h), axis=-1)
    
    for i in range(3):
        # d_j Atilde^ij
        divA = dAtilde[..., i, 0, 0] + dAtilde[..., i, 1, 1] + dAtilde[..., i, 2, 2]
        # 6 Atilde^ij d_j phi
        term2 = 6.0 * (Atilde[..., i, 0] * dphi[..., 0] +
                       Atilde[..., i, 1] * dphi[..., 1] +
                       Atilde[..., i, 2] * dphi[..., 2])
        # -2/3 d^i K
        term3 = -2.0/3.0 * (ginv[..., i, 0] * dK[..., 0] +
                            ginv[..., i, 1] * dK[..., 1] +
                            ginv[..., i, 2] * dK[..., 2])
        M[..., i] = divA + term2 + term3 - 8.0 * np.pi * S_i[..., i]
    
    return M


__all__ = [
    "apply_dirichlet_bc",
    "apply_neumann_bc",
    "apply_sommerfeld_bc",
    "apply_constraint_preserving_bc_scalar",
    "apply_constraint_preserving_bc_complex",
    "apply_constraint_preserving_bc_ptum",
    "hamiltonian_constraint_residual",
    "momentum_constraint_residual",
]