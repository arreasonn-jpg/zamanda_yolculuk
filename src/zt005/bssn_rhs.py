"""BSSN RHS: conformal Christoffel, Ricci, constraint damping.

Operates pointwise. Derivatives are caller-provided (via FD on a grid
or analytic). This keeps the module grid-agnostic and testable.
"""
from __future__ import annotations
import numpy as np
from .bssn import BSSN, euler_step


def conformal_christoffel(gtilde, dgtilde):
    """Gamma^i_jk = 1/2 gt^{il} (d_j gt_lk + d_k gt_lj - d_l gt_jk).
    gtilde: (3,3). dgtilde: (3,3,3) with dgtilde[i,j,k] = d_k gt_ij."""
    gt = np.asarray(gtilde, float)
    dg = np.asarray(dgtilde, float)
    gi = np.linalg.inv(gt)
    G = np.zeros((3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                s = 0.0
                for l in range(3):
                    s += gi[i, l] * (dg[j, l, k] + dg[k, l, j] - dg[l, j, k])
                G[i, j, k] = 0.5 * s
    return G


def conformal_connection(G):
    """Gamma^i = gt^{jk} Gamma^i_jk -> shape (3,)."""
    return np.einsum('ijk,jk->i', G, np.eye(3))


def ricci_tensor(G, dG):
    """R_ij = d_k Gamma^k_ij - d_j Gamma^k_ik
              + Gamma^k_kl Gamma^l_ij - Gamma^k_jl Gamma^l_ik.
    dG[k,i,j,l] = d_l Gamma^k_ij."""
    R = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            s = 0.0
            for k in range(3):
                s += dG[k, i, j, k] - dG[k, i, k, j]
                for l in range(3):
                    s += G[k, k, l] * G[l, i, j] - G[k, j, l] * G[l, i, k]
            R[i, j] = s
    return R


def constraint_damping_term(H, M, eta=0.1):
    """Adds -eta H to dK and -eta M^i to dGamma^i (Z4c-inspired, minimal)."""
    dK_damp = -eta * float(H)
    dGamma_damp = -eta * np.asarray(M, float)
    return dK_damp, dGamma_damp


def bssn_rhs(b: BSSN, rho: float, S: np.ndarray,
             dgtilde: np.ndarray, dG: np.ndarray,
             eta: float = 0.1) -> dict:
    """Assemble RHS dict consumable by bssn.euler_step.
    dgtilde[i,j,k] = d_k gtilde_ij ; dG[k,i,j,l] = d_l Gamma^k_ij."""
    from .bssn import hamiltonian_constraint, momentum_constraint
    G = conformal_christoffel(b.gtilde, dgtilde)
    Rij = ricci_tensor(G, dG)
    H = hamiltonian_constraint(b, rho)
    M = momentum_constraint(b, S)
    dK_damp, dGamma_damp = constraint_damping_term(H, M, eta)
    # Skeleton RHS: evolution of Atilde ~ -2 alpha R^TF ; phi ~ -1/6 alpha K
    alpha = 1.0
    Atilde_dot = -2.0 * alpha * (Rij - (np.trace(Rij) / 3.0) * b.gtilde)
    return {
        "dphi": -alpha * b.K / 6.0,
        "dgt": np.zeros_like(b.gtilde),
        "dK": -alpha * float(np.sum(Rij)) + dK_damp,
        "dAt": Atilde_dot,
        "dGamma": dGamma_damp,
    }


__all__ = ["conformal_christoffel", "conformal_connection",
           "ricci_tensor", "constraint_damping_term", "bssn_rhs"]
