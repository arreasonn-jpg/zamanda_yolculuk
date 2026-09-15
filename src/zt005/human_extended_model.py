"""Human Extended-Body Proxy Model & Tidal-Force Safety Gate.

Addresses Section 16 and 17:
- Human extended-body proxy:
    mass = 75 kg, height = 1.75 m, width = 0.45 m, depth = 0.25 m
    center-of-mass worldline x^mu(tau), 4-velocity u^mu, orientation triad
- Tidal field & geodesic deviation:
    Electric part of Riemann tensor in traveler rest frame:
        E_ij = R_{0_hat i 0_hat j} = R^mu_anb e^a_(0) e^b_(0) e_(i)mu e_(j)^nu
    Tidal acceleration:
        a^i_tidal = - E_ij xi^j
    Jerk:
        j^i = d a^i_tidal / d tau
- Safety Gate:
    Pass if head-to-toe tidal acceleration <= 10 g (98.1 m/s^2)
    Fail if tidal force causes structural injury / spaghettification.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

C = 299_792_458.0
G_EARTH = 9.80665  # 1 g [m/s^2]


@dataclass(frozen=True)
class HumanProxy:
    mass_kg: float = 75.0
    height_m: float = 1.75   # Head-to-toe separation
    width_m: float = 0.45    # Shoulder-to-shoulder separation
    depth_m: float = 0.25    # Chest-to-back separation


def construct_orthonormal_tetrad(g: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Construct an orthonormal tetrad e_(a)^mu with e_(0)^mu = u^mu, satisfying g(e_a, e_b) = eta_ab."""
    g = np.asarray(g, float)
    u = np.asarray(u, float)
    norm = float(u @ g @ u)
    if norm >= 0:
        raise ValueError(f"u must be timelike, got norm = {norm}")

    e0 = u / np.sqrt(abs(norm))

    # Gram-Schmidt orthogonalization for spatial triad e1, e2, e3
    tetrad = [e0]
    candidates = np.eye(4)[1:]  # (0,1,0,0), (0,0,1,0), (0,0,0,1)

    for cand in candidates:
        v = cand.copy()
        for e in tetrad:
            eta_val = -1.0 if len(tetrad) == 1 else 1.0
            proj = float(v @ g @ e)
            v -= proj / eta_val * e
        v_norm = float(v @ g @ v)
        if abs(v_norm) < 1e-12:
            # Degenerate choice, perturb
            v = cand + np.array([0.0, 0.1, 0.1, 0.1])
            for e in tetrad:
                eta_val = -1.0 if len(tetrad) == 1 else 1.0
                proj = float(v @ g @ e)
                v -= proj / eta_val * e
            v_norm = float(v @ g @ v)
        v = v / np.sqrt(abs(v_norm))
        tetrad.append(v)

    return np.asarray(tetrad)  # shape (4, 4), tetrad[a, mu] = e_(a)^mu


def compute_tidal_tensor(riemann: np.ndarray, tetrad: np.ndarray) -> np.ndarray:
    """Compute the electric tidal tensor E_ij = R_{0_hat i 0_hat j} in traveler's rest frame.

    E_ij = R^mu_{alpha nu beta} e_(i)_mu e_(0)^alpha e_(j)^nu e_(0)^beta
    """
    R = np.asarray(riemann, float)
    e = np.asarray(tetrad, float)

    # e[0] is e_(0), e[1..3] are e_(1..3)
    e0 = e[0]
    E = np.zeros((3, 3), float)

    for i in range(3):
        ei = e[i + 1]
        for j in range(3):
            ej = e[j + 1]
            # R^m_anb ei_m e0^a ej^n e0^b
            val = np.einsum("manb,m,a,n,b->", R, ei, e0, ej, e0)
            E[i, j] = val

    # Symmetrize
    return 0.5 * (E + E.T)


def evaluate_human_safety_gate(
    riemann: np.ndarray,
    g: np.ndarray,
    u: np.ndarray,
    human: HumanProxy | None = None,
    max_accel_g: float = 10.0,
    max_jerk_m_s3: float = 50.0,
    dE_dtau: np.ndarray | None = None,
) -> dict:
    """Evaluate tidal forces on human body proxy and determine Safety Gate PASS/FAIL."""
    human = human or HumanProxy()
    tetrad = construct_orthonormal_tetrad(g, u)
    E_ij = compute_tidal_tensor(riemann, tetrad)

    # Head-to-toe separation vector along e_(3) (height axis)
    xi_head_toe = np.array([0.0, 0.0, human.height_m])
    # Tidal acceleration vector: a^i = - E_ij xi^j [in geometric units c=1, SI: c^2 * E_ij * xi^j]
    # Note: Riemann in SI has units [1/m^2], so c^2 * R has units [1/s^2], and a_tidal = c^2 * E * xi [m/s^2]
    a_tidal = (C**2) * (E_ij @ xi_head_toe)
    a_tidal_mag = float(np.linalg.norm(a_tidal))
    a_tidal_in_g = float(a_tidal_mag / G_EARTH)

    # Chest-to-back tidal compression
    xi_depth = np.array([human.depth_m, 0.0, 0.0])
    a_depth = (C**2) * (E_ij @ xi_depth)
    a_depth_mag = float(np.linalg.norm(a_depth))

    # Jerk
    if dE_dtau is not None:
        jerk = (C**2) * (dE_dtau @ xi_head_toe)
        jerk_mag = float(np.linalg.norm(jerk))
    else:
        jerk_mag = 0.0

    # Gate logic
    passed_accel = bool(a_tidal_in_g <= max_accel_g)
    passed_jerk = bool(jerk_mag <= max_jerk_m_s3)
    safety_pass = bool(passed_accel and passed_jerk)

    if safety_pass:
        verdict = "PASS: Tidal forces are structurally survivable for a human subject."
    else:
        verdict = f"FAIL: Lethal tidal forces detected ({a_tidal_in_g:.2e} g > {max_accel_g} g threshold) — subject would suffer structural spaghettification."

    return {
        "safety_gate_pass": safety_pass,
        "verdict": verdict,
        "head_to_toe_tidal_acceleration_m_s2": a_tidal_mag,
        "head_to_toe_tidal_acceleration_g": a_tidal_in_g,
        "max_allowable_g": max_accel_g,
        "chest_compression_acceleration_m_s2": a_depth_mag,
        "jerk_m_s3": jerk_mag,
        "tidal_tensor_E_ij": E_ij.tolist(),
        "human_parameters": {
            "mass_kg": human.mass_kg,
            "height_m": human.height_m,
            "width_m": human.width_m,
            "depth_m": human.depth_m,
        },
    }
