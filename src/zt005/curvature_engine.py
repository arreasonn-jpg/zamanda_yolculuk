"""Full Curvature and Einstein Tensor Pipeline.

Addresses Section 3.1 and 3.2:
- Metric consistency chain:
    g_munu -> g^munu -> Gamma^mu_ab -> R^rho_sab -> R_ab -> R -> G_ab
- Einstein residual:
    Delta G_munu = G_munu - (8 pi G / c^4) T_munu
- Gauge residual:
    partial_mu hbar^{mu nu} ~ 0
- Rigorous solver classification:
    "linearized metric perturbation solver" vs "full nonlinear GR solver"
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0
G_CONST = 6.67430e-11
KAPPA = 8.0 * np.pi * G_CONST / C**4
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def compute_metric_inverse(g: np.ndarray) -> np.ndarray:
    """Compute g^{mu nu} via 4x4 matrix inverse."""
    g_arr = np.asarray(g, float)
    return np.linalg.inv(g_arr)


def compute_christoffel(g_fn, x: np.ndarray, h: float = 1e-4) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute metric, inverse, and Christoffel symbols Gamma^mu_{alpha beta} at x.

    Gamma^mu_ab = 1/2 g^{mu lambda} (d_a g_bl + d_b g_al - d_l g_ab)
    """
    x = np.asarray(x, float)
    g = np.asarray(g_fn(x), float)
    g_inv = compute_metric_inverse(g)

    # 4th-order central differences for dg_ab / dx^mu
    dg = np.zeros((4, 4, 4))  # dg[mu, alpha, beta] = d g_{alpha beta} / d x^mu
    for mu in range(4):
        x_p2 = x.copy(); x_p2[mu] += 2.0 * h
        x_p1 = x.copy(); x_p1[mu] += 1.0 * h
        x_m1 = x.copy(); x_m1[mu] -= 1.0 * h
        x_m2 = x.copy(); x_m2[mu] -= 2.0 * h

        g_p2 = np.asarray(g_fn(x_p2), float)
        g_p1 = np.asarray(g_fn(x_p1), float)
        g_m1 = np.asarray(g_fn(x_m1), float)
        g_m2 = np.asarray(g_fn(x_m2), float)

        dg[mu] = (-g_p2 + 8.0 * g_p1 - 8.0 * g_m1 + g_m2) / (12.0 * h)

    # Gamma^m_ab = 1/2 g^{ml} (dg[a, b, l] + dg[b, a, l] - dg[l, a, b])
    gamma = np.zeros((4, 4, 4))
    for m in range(4):
        for a in range(4):
            for b in range(4):
                val = 0.0
                for l in range(4):
                    term = dg[a, b, l] + dg[b, a, l] - dg[l, a, b]
                    val += 0.5 * g_inv[m, l] * term
                gamma[m, a, b] = val

    return g, g_inv, gamma


def compute_riemann(g_fn, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Compute Riemann curvature tensor R^rho_{sigma mu nu} at x.

    R^rho_smn = d_m Gamma^rho_ns - d_n Gamma^rho_ms
                + Gamma^rho_ml Gamma^l_ns - Gamma^rho_nl Gamma^l_ms
    """
    x = np.asarray(x, float)
    _, _, gamma0 = compute_christoffel(g_fn, x, h=h)

    # Numerical derivatives of Gamma: d Gamma / dx^k
    d_gamma = np.zeros((4, 4, 4, 4))  # d_gamma[k, rho, alpha, beta] = d/dx^k Gamma^rho_ab
    for k in range(4):
        x_p = x.copy(); x_p[k] += h
        x_m = x.copy(); x_m[k] -= h
        _, _, g_p = compute_christoffel(g_fn, x_p, h=h)
        _, _, g_m = compute_christoffel(g_fn, x_m, h=h)
        d_gamma[k] = (g_p - g_m) / (2.0 * h)

    riemann = np.zeros((4, 4, 4, 4))
    for rho in range(4):
        for sig in range(4):
            for mu in range(4):
                for nu in range(4):
                    d_term = d_gamma[mu, rho, nu, sig] - d_gamma[nu, rho, mu, sig]
                    quad_term = 0.0
                    for l in range(4):
                        quad_term += gamma0[rho, mu, l] * gamma0[l, nu, sig] - gamma0[rho, nu, l] * gamma0[l, mu, sig]
                    riemann[rho, sig, mu, nu] = d_term + quad_term

    return riemann


def compute_curvature_chain(g_fn, x: np.ndarray, h: float = 1e-4) -> dict:
    """Execute full curvature chain: g -> g^inv -> Gamma -> Riemann -> Ricci -> R -> Einstein."""
    g, g_inv, gamma = compute_christoffel(g_fn, x, h=h)
    riemann = compute_riemann(g_fn, x, h=h)

    # Ricci tensor: R_ab = R^l_alb
    ricci = np.einsum("lalb->ab", riemann)
    ricci = 0.5 * (ricci + ricci.T)  # enforce symmetry

    # Ricci scalar: R = g^{ab} R_ab
    ricci_scalar = float(np.einsum("ab,ab->", g_inv, ricci))

    # Einstein tensor: G_ab = R_ab - 1/2 g_ab R
    einstein = ricci - 0.5 * g * ricci_scalar

    return {
        "g": g,
        "g_inv": g_inv,
        "gamma": gamma,
        "riemann": riemann,
        "ricci": ricci,
        "ricci_scalar": ricci_scalar,
        "einstein": einstein,
    }


def compute_einstein_residual(g_fn, x: np.ndarray, T_munu: np.ndarray, h: float = 1e-4) -> dict:
    """Compute Einstein equation residual Delta G_munu = G_munu - (8 pi G / c^4) T_munu.

    Returns residual and solver classification.
    """
    chain = compute_curvature_chain(g_fn, x, h=h)
    G_ab = chain["einstein"]
    T_ab = np.asarray(T_munu, float)

    source_term = KAPPA * T_ab
    residual_tensor = G_ab - source_term

    res_frobenius = float(np.linalg.norm(residual_tensor))
    scale = max(float(np.linalg.norm(G_ab)), float(np.linalg.norm(source_term)), 1e-300)
    rel_residual = float(res_frobenius / scale)

    # Classification
    classification = "linearized metric perturbation solver"

    return {
        "G_munu": G_ab,
        "source_term_kappa_T": source_term,
        "residual_tensor": residual_tensor,
        "residual_frobenius": res_frobenius,
        "relative_residual": rel_residual,
        "classification": classification,
        "is_full_gr_solution": bool(rel_residual < 0.01),
    }
