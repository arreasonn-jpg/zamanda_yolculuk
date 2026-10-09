"""BSSN-form 3+1 GR skeleton (minimal, constraint-checking only).

Not a full numerical-relativity solver. Provides:
- ADM 3+1 variable container
- BSSN conformal variable container
- Hamiltonian & momentum constraint residuals
- One explicit Euler step of the BSSN evolution system
- Linear-regime sanity vs linearized_gr
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

C = 2.99792458e8
G = 6.67430e-11


@dataclass
class ADM:
    alpha: float
    beta: np.ndarray      # shape (3,)
    gamma: np.ndarray     # shape (3,3)
    K: np.ndarray         # shape (3,3)


@dataclass
class BSSN:
    phi: float            # conformal factor: gamma_ij = e^{4 phi} gtilde_ij
    gtilde: np.ndarray    # det = 1
    K: float              # trace of extrinsic curvature
    Atilde: np.ndarray    # trace-free part
    Gamma: np.ndarray     # conformal connection, shape (3,)


def adm_to_bssn(adm: ADM) -> BSSN:
    g = np.asarray(adm.gamma, float)
    detg = np.linalg.det(g)
    phi = (1.0 / 12.0) * np.log(max(detg, 1e-300))
    gt = np.exp(-4.0 * phi) * g
    trK = float(np.trace(adm.K @ np.linalg.inv(g)))
    A = adm.K - (trK / 3.0) * g
    At = np.exp(-4.0 * phi) * A
    return BSSN(phi=float(phi), gtilde=gt, K=trK, Atilde=At,
                Gamma=np.zeros(3))


def bssn_to_adm(b: BSSN, alpha: float, beta: np.ndarray) -> ADM:
    g = np.exp(4.0 * b.phi) * b.gtilde
    K = np.exp(4.0 * b.phi) * b.Atilde + (b.K / 3.0) * g
    return ADM(alpha=alpha, beta=np.asarray(beta, float), gamma=g, K=K)


def hamiltonian_constraint(b: BSSN, rho: float) -> float:
    """H = R + K^2 - K_ij K^ij - 16 pi G rho / c^2  (simplified R)."""
    R = 0.0  # flat at this skeleton level
    K2 = b.K ** 2
    KijKij = float(np.sum(b.Atilde ** 2))
    return float(R + K2 - KijKij - 16.0 * np.pi * G * rho / C ** 2)


def momentum_constraint(b: BSSN, S: np.ndarray) -> np.ndarray:
    """M^i = D_j (K^ij - gamma^ij K) - 8 pi G S^i / c^4 (skeleton: 0)."""
    return -8.0 * np.pi * G * np.asarray(S, float) / C ** 4


def euler_step(b: BSSN, rhs: dict, dt: float) -> BSSN:
    """One explicit Euler update with caller-provided RHS dict."""
    return BSSN(
        phi=b.phi + dt * rhs.get("dphi", 0.0),
        gtilde=b.gtilde + dt * rhs.get("dgt", np.zeros_like(b.gtilde)),
        K=b.K + dt * rhs.get("dK", 0.0),
        Atilde=b.Atilde + dt * rhs.get("dAt", np.zeros_like(b.Atilde)),
        Gamma=b.Gamma + dt * rhs.get("dGamma", np.zeros_like(b.Gamma)),
    )


def linear_regime_check(h: float, tol: float = 1e-3) -> bool:
    """True if h is small enough that BSSN ~ linearized GR."""
    return abs(h) < tol


__all__ = ["ADM", "BSSN", "adm_to_bssn", "bssn_to_adm",
           "hamiltonian_constraint", "momentum_constraint",
           "euler_step", "linear_regime_check"]
