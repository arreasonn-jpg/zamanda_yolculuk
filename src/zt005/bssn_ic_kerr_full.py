"""Full Kerr IC via Kerr-Schild Cartesian coordinates (v3.0-alpha1).

KS form:
  ds^2 = -(1-2H) dt^2 + 4H k_i dt dx^i + (delta_ij + 2H k_i k_j) dx^i dx^j
  H = M r^3 / (r^4 + a^2 z^2)
  k_i = ((r x + a y)/(r^2+a^2), (r y - a x)/(r^2+a^2), z/r)
  r^2 = (rho^2 - a^2 + sqrt((rho^2-a^2)^2 + 4a^2 z^2))/2

Works for |a| <= M. Horizon-penetrating.

3+1:
  alpha = 1/sqrt(1 + 2H)
  beta^i = 2H k^i / (1 + 2H)
  gamma_ij = delta_ij + 2H k_i k_j
  K_ij = (1/(2 alpha)) (D_i beta_j + D_j beta_i)
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_grid import d_central


def kerr_schild_metric(N, L, M=1.0, a=0.0):
    """Return (r, H, k_i, alpha, beta_i, beta^i, gamma_ij) on (N,N,N) grid."""
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho2 = X**2 + Y**2 + Z**2
    disc = np.sqrt(np.maximum((rho2 - a**2)**2 + 4.0*a**2*Z**2, 0.0))
    r2 = 0.5 * (rho2 - a**2 + disc)
    r2 = np.maximum(r2, 1e-12)
    R = np.sqrt(r2)
    H = M * R**3 / np.maximum(R**4 + a**2 * Z**2, 1e-30)
    den = np.maximum(R**2 + a**2, 1e-30)
    kx = (R * X + a * Y) / den
    ky = (R * Y - a * X) / den
    kz = Z / np.maximum(R, 1e-30)
    k = np.stack([kx, ky, kz], axis=-1)
    alpha = 1.0 / np.sqrt(1.0 + 2.0 * H)
    beta_i = 2.0 * H[..., None] * k           # covariant shift
    # gamma_ij = delta_ij + 2H k_i k_j
    gamma = np.broadcast_to(np.eye(3)[None,None,None,:,:],
                            (N,N,N,3,3)).copy()
    gamma = gamma + 2.0 * H[..., None, None] * np.einsum('...i,...j->...ij', k, k)
    ginv = np.linalg.inv(gamma)
    beta_up = np.einsum('...ij,...j->...i', ginv, beta_i)
    return R, H, k, alpha, beta_i, beta_up, gamma


def kerr_schild_state(N, L, M=1.0, a=0.0):
    """Build BSSNState for Kerr (KS Cartesian). |a| <= M."""
    if abs(a) > M:
        raise ValueError("|a| <= M required")
    R, H, k, alpha, beta_i, beta_up, gamma = kerr_schild_metric(N, L, M, a)
    h = L / (N - 1)
    # BSSN conformal decomposition
    detg = np.linalg.det(gamma)
    phi = (1.0/12.0) * np.log(np.maximum(detg, 1e-30))
    gtilde = np.exp(-4.0 * phi[..., None, None]) * gamma
    ginv = np.linalg.inv(gamma)
    # Extrinsic curvature via FD: K_ij = (1/(2 a)) (D_i beta_j + D_j beta_i)
    dbeta = np.zeros((N, N, N, 3, 3))     # dbeta[...,i,j] = d_j beta_i
    for i in range(3):
        comp = beta_i[..., i]
        dx, dy, dz = d_central(comp, h)
        dbeta[..., i, 0] = dx
        dbeta[..., i, 1] = dy
        dbeta[..., i, 2] = dz
    dgamma = np.zeros((N, N, N, 3, 3, 3))  # dgamma[...,i,j,k] = d_k gamma_ij
    for i in range(3):
        for j in range(3):
            comp = gamma[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dgamma[..., i, j, 0] = dx
            dgamma[..., i, j, 1] = dy
            dgamma[..., i, j, 2] = dz
    # Christoffel
    Gamma = np.zeros((N, N, N, 3, 3, 3))
    for kk in range(3):
        for i in range(3):
            for j in range(3):
                s = np.zeros((N, N, N))
                for l in range(3):
                    s += 0.5 * ginv[..., kk, l] * (
                        dgamma[..., l, j, i] + dgamma[..., l, i, j]
                        - dgamma[..., i, j, l])
                Gamma[..., kk, i, j] = s
    # D_i beta_j
    Dbeta = dbeta.copy()
    for i in range(3):
        for j in range(3):
            for kk in range(3):
                Dbeta[..., i, j] -= Gamma[..., kk, i, j] * beta_i[..., kk]
    K_ij = np.zeros((N, N, N, 3, 3))
    for i in range(3):
        for j in range(3):
            K_ij[..., i, j] = (Dbeta[..., i, j] + Dbeta[..., j, i]) / (2.0 * alpha)
    K_tr = np.einsum('...ij,...ij->...', ginv, K_ij)
    Atilde = np.exp(-4.0 * phi[..., None, None]) * (
        K_ij - (1.0/3.0) * gamma * K_tr[..., None, None])
    return BSSNState(
        phi=phi, gtilde=gtilde, K=K_tr, Atilde=Atilde,
        Gamma=np.zeros((N, N, N, 3)),
        alpha=alpha, beta=beta_up,
    )


__all__ = ["kerr_schild_metric", "kerr_schild_state"]
