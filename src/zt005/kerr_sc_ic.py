"""Kerr + SC smoothed core IC (v6.3-alpha1).

Combines Kerr-Schild Cartesian metric with a smoothed mass profile:
  M_eff(r) = M (1 - e^{-r/rc} (1 + r/rc))
so that M_eff(0) = 0 (regular center) and M_eff(r >> rc) = M
(recover Kerr).

This produces a rotating regular-core metric: SC-like regularization
of the Kerr singularity at r=0, while preserving the ring structure
at larger r (Kerr-Schild coordinates).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_grid import d_central


def kerr_sc_metric(N, L, M=1.0, a=0.0, rc=1.0):
    """Kerr-Schild metric with smoothed mass M_eff(r).

    Returns dict with gamma_ij, alpha, beta_i (covariant), beta_up.
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho2 = X**2 + Y**2 + Z**2
    # Smoothed mass
    R = np.sqrt(rho2) + 1e-30
    M_eff = M * (1.0 - np.exp(-R/rc) * (1.0 + R/rc))
    # Kerr-Schild functions with M_eff
    disc = np.sqrt(np.maximum((rho2 - a**2)**2 + 4.0*a**2*Z**2, 0.0))
    r2 = 0.5 * (rho2 - a**2 + disc)
    r2 = np.maximum(r2, 1e-12)
    Rks = np.sqrt(r2)
    H = M_eff * Rks**3 / np.maximum(Rks**4 + a**2 * Z**2, 1e-30)
    den = np.maximum(Rks**2 + a**2, 1e-30)
    kx = (Rks * X + a * Y) / den
    ky = (Rks * Y - a * X) / den
    kz = Z / np.maximum(Rks, 1e-30)
    k = np.stack([kx, ky, kz], axis=-1)
    alpha = 1.0 / np.sqrt(1.0 + 2.0 * H)
    beta_i = 2.0 * H[..., None] * k
    gamma = np.broadcast_to(np.eye(3)[None,None,None,:,:],
                            (N,N,N,3,3)).copy()
    gamma = gamma + 2.0 * H[..., None, None] * np.einsum('...i,...j->...ij', k, k)
    ginv = np.linalg.inv(gamma)
    beta_up = np.einsum('...ij,...j->...i', ginv, beta_i)
    return {
        "gamma": gamma, "alpha": alpha, "beta_i": beta_i, "beta_up": beta_up,
        "H": H, "M_eff": M_eff, "R": R,
    }


def kerr_sc_state(N, L, M=1.0, a=0.0, rc=1.0):
    """Build BSSNState for Kerr + SC smoothed core."""
    if abs(a) > M:
        raise ValueError("|a| <= M required")
    met = kerr_sc_metric(N, L, M=M, a=a, rc=rc)
    gamma = met["gamma"]
    alpha = met["alpha"]
    beta_i = met["beta_i"]
    beta_up = met["beta_up"]
    h = L / (N - 1)
    detg = np.linalg.det(gamma)
    phi = (1.0/12.0) * np.log(np.maximum(detg, 1e-30))
    gtilde = np.exp(-4.0 * phi[..., None, None]) * gamma
    ginv = np.linalg.inv(gamma)
    # Extrinsic curvature: K_ij = (1/(2 a)) (D_i beta_j + D_j beta_i)
    dbeta = np.zeros((N, N, N, 3, 3))
    for i in range(3):
        comp = beta_i[..., i]
        dx, dy, dz = d_central(comp, h)
        dbeta[..., i, 0] = dx
        dbeta[..., i, 1] = dy
        dbeta[..., i, 2] = dz
    dgamma = np.zeros((N, N, N, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            comp = gamma[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dgamma[..., i, j, 0] = dx
            dgamma[..., i, j, 1] = dy
            dgamma[..., i, j, 2] = dz
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


__all__ = ["kerr_sc_metric", "kerr_sc_state"]
