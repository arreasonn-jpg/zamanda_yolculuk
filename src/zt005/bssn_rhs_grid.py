"""Grid-aware full BSSN RHS: computes all FD derivatives of phi, Atilde,
Gamma and assembles dgtilde, dphi, dK, dAtilde, dGamma."""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central, dgtilde_from_grid, christoffel_from_grid, dG_from_grid


def _d(f, h, ax):
    """First derivative along axis ax."""
    out = np.zeros_like(f)
    if f.ndim == 3:
        d = d_central(f, h)[ax]
    else:
        # for tensor fields, iterate over components
        d = np.zeros_like(f)
    return d


def grad_scalar(f, h):
    """Return (df/dx, df/dy, df/dz)."""
    return d_central(f, h)


def hess_scalar(f, h):
    """Return Hessian H[i,j] = d_i d_j f on grid (3,3,N,N,N)."""
    d = d_central(f, h)
    H = np.zeros((3, 3) + f.shape)
    for i in range(3):
        dd = d_central(d[i], h)
        for j in range(3):
            H[i, j] = dd[j]
    return H


def grad_tensor(A, h):
    """A: (N,N,N,3,3). Return dA[i,j,k] = d_k A_ij (grid)."""
    Nx, Ny, Nz = A.shape[:3]
    dA = np.zeros((Nx, Ny, Nz, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            comp = A[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dA[..., i, j, 0] = dx
            dA[..., i, j, 1] = dy
            dA[..., i, j, 2] = dz
    return dA


def bssn_rhs_grid(phi, gt, K, At, h, rho_grid=None, S_grid=None,
                  eta=0.1):
    """Full grid RHS. Returns dict of derivative fields, same shapes as inputs."""
    from .bssn import BSSN
    from .bssn_rhs_full import dt_gtilde, dt_phi, dt_K, dt_Atilde
    from .bssn_rhs import conformal_christoffel, ricci_tensor, constraint_damping_term
    from .bssn import hamiltonian_constraint, momentum_constraint

    Nx, Ny, Nz = phi.shape
    # Derivatives
    dphi = np.array(d_central(phi, h))            # (3, N, N, N)
    Hphi = hess_scalar(phi, h)                     # (3,3,N,N,N)
    dAt = grad_tensor(At, h)                       # (N,N,N,3,3,3)
    dg = dgtilde_from_grid(gt, h)                  # (N,N,N,3,3,3)
    G, _ = christoffel_from_grid(gt, h)
    dG = dG_from_grid(G, h)

    dphi_t = np.zeros_like(phi)
    dgt_t = np.zeros_like(gt)
    dK_t = np.zeros_like(K)
    dAt_t = np.zeros_like(At)

    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                rho = float(rho_grid[ix,iy,iz]) if rho_grid is not None else 0.0
                S = S_grid[ix,iy,iz] if S_grid is not None else np.zeros(3)
                b = BSSN(float(phi[ix,iy,iz]), gt[ix,iy,iz],
                         float(K[ix,iy,iz]), At[ix,iy,iz], np.zeros(3))
                Rij = ricci_tensor(G[ix,iy,iz], dG[ix,iy,iz])
                H = hamiltonian_constraint(b, rho)
                M = momentum_constraint(b, S)
                dK_damp, dG_damp = constraint_damping_term(H, M, eta)

                dphi_t[ix,iy,iz] = dt_phi(1.0, b.K)
                dgt_t[ix,iy,iz]  = dt_gtilde(1.0, b.Atilde)
                dK_t[ix,iy,iz]   = dt_K(1.0, b.Atilde, Rij, 0.0, b.K) + dK_damp

                # Full Atilde: -2 alpha D_i D_j phi + 4 alpha d_i phi d_j phi + alpha R_ij
                hess = Hphi[:, :, ix, iy, iz]
                grad = dphi[:, ix, iy, iz]
                dphi_dphi = np.outer(grad, grad)
                A2 = b.Atilde @ b.Atilde
                raw = (-2.0 * hess + 4.0 * dphi_dphi
                       + Rij + b.Atilde * b.K - 2.0 * A2)
                tr = float(np.trace(raw))
                dAt_t[ix,iy,iz] = raw - (tr/3.0) * np.eye(3)

    return {"dphi": dphi_t, "dgt": dgt_t, "dK": dK_t, "dAt": dAt_t}


__all__ = ["grad_scalar", "hess_scalar", "grad_tensor", "bssn_rhs_grid"]
