"""Full BSSN RHS from ADM (partial: dgtilde, dphi, dK, dAt, dGamma).

Based on standard BSSN equations (Alcubierre 2008, Sec. 11.3).
Assumes: alpha=1 (lapse), beta=0 (shift) for this skeleton; these will
be promoted to fields in a later stage.
"""
from __future__ import annotations
import numpy as np
from .bssn import BSSN
from .bssn_rhs import conformal_christoffel, ricci_tensor, constraint_damping_term
from .bssn import hamiltonian_constraint, momentum_constraint


def dt_gtilde(alpha, Atilde):
    """d_t gtilde_ij = -2 alpha Atilde_ij  (+ L_beta gtilde = 0)."""
    return -2.0 * alpha * np.asarray(Atilde, float)


def dt_phi(alpha, K):
    """d_t phi = -alpha K / 6."""
    return -alpha * K / 6.0


def dt_K(alpha, Atilde, R, D2alpha, K):
    """d_t K = -D^2 alpha + alpha (Atilde_ij Atilde^ij + K^2/3) + 4 pi alpha (rho + S).
    Here we skip the source term and D^2 alpha (alpha=1)."""
    A2 = float(np.sum(Atilde * Atilde))
    return alpha * (A2 + K ** 2 / 3.0)


def dt_Atilde(alpha, phi, Atilde, R, K):
    """d_t Atilde_ij = e^{-4phi} [ -2 alpha D_i D_j phi + 4 alpha d_i phi d_j phi
        + alpha R_ij - 2 alpha Atilde_ik Atilde^k_j
        + alpha Atilde_ij K ]^TF  (trace-free part)."""
    A = np.asarray(Atilde, float)
    Rij = np.asarray(R, float)
    A2 = A @ A
    full = alpha * (Rij + A * K) - 2.0 * alpha * A2
    tr = float(np.trace(full))
    return full - (tr / 3.0) * np.eye(3)


def dt_Gamma(alpha, Atilde, dAtilde, dgtilde, phi):
    """d_t Gamma^i ~ -2 d_j (alpha Atilde^ij) + 2 alpha Atilde^ij d_j phi
    (skeleton; full version has gamma-tilde-trace terms)."""
    # dAtilde[j,i,j] contraction: d_j Atilde^ji
    dA = np.asarray(dAtilde, float)  # shape (3,3,3): dA[i,j,k]=d_k Atilde_ij
    # We use Atilde @ ... approximate here; full form in P10b.
    return -2.0 * np.einsum('jij->i', dA)


def bssn_rhs_full(b: BSSN, rho, S, dgtilde, dG, eta=0.1) -> dict:
    """Assemble full RHS dict (drop-in for bssn.euler_step)."""
    alpha = 1.0
    G = conformal_christoffel(b.gtilde, dgtilde)
    R = ricci_tensor(G, dG)
    H = hamiltonian_constraint(b, rho)
    M = momentum_constraint(b, S)
    dK_damp, dGamma_damp = constraint_damping_term(H, M, eta)
    dgt = dt_gtilde(alpha, b.Atilde)
    dphi = dt_phi(alpha, b.K)
    dK = dt_K(alpha, b.Atilde, R, 0.0, b.K) + dK_damp
    # dAtilde needs spatial derivatives of Atilde (approximated by 0 here)
    dAt = dt_Atilde(alpha, b.phi, b.Atilde, R, b.K)
    dGamma = np.zeros(3) + dGamma_damp
    return {"dphi": dphi, "dgt": dgt, "dK": dK, "dAt": dAt, "dGamma": dGamma}


__all__ = ["dt_gtilde", "dt_phi", "dt_K", "dt_Atilde", "dt_Gamma",
           "bssn_rhs_full"]
