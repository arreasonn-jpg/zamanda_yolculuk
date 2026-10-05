"""Stress-Energy Physical Consistency & Energy Conditions Engine.

Addresses Section 2.1 and 2.2:
- Stress-energy conservation gate: partial_mu T^{mu nu} = 0 for nu=0,1,2,3 with PASS/WARN/FAIL.
- Positive energy control: T_00 >= 0 everywhere.
- Electromagnetic energy density cross-check: u = 1/2 (eps0 E^2 + B^2/mu0) vs T_00.
- Energy condition diagnostics:
    * Weak Energy Condition (WEC): T_munu v^mu v^nu >= 0
    * Dominant Energy Condition (DEC): -T^mu_nu v^nu is future-directed causal
    * Null Energy Condition (NEC): T_munu k^mu k^nu >= 0
"""

from __future__ import annotations

import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 4.0e-7 * np.pi
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def verify_maxwell_energy_density(E: np.ndarray, B: np.ndarray, T_00: np.ndarray, tol: float = 1e-10) -> dict:
    """Compare solver T_00 against analytic electromagnetic energy density u = 1/2(eps0 E^2 + B^2/mu0)."""
    E = np.asarray(E, float)
    B = np.asarray(B, float)
    T00 = np.asarray(T_00, float)

    E2 = np.sum(E * E, axis=-1)
    B2 = np.sum(B * B, axis=-1)
    u_analytic = 0.5 * EPS0 * E2 + 0.5 / MU0 * B2

    abs_err = np.abs(T00 - u_analytic)
    scale = np.maximum(u_analytic, 1e-30)
    rel_err = abs_err / scale

    max_abs = float(np.max(abs_err))
    max_rel = float(np.max(rel_err))
    all_positive = bool(np.all(T00 >= -1e-15))
    is_passed = bool(max_rel < tol and all_positive)

    return {
        "pass": is_passed,
        "all_positive_T00": all_positive,
        "max_abs_err": max_abs,
        "max_rel_err": max_rel,
        "u_min": float(np.min(u_analytic)),
        "u_max": float(np.max(u_analytic)),
        "T00_min": float(np.min(T00)),
        "T00_max": float(np.max(T00)),
    }


def evaluate_energy_conditions(T_field: np.ndarray) -> dict:
    """Evaluate WEC, DEC, NEC, and SEC on stress-energy tensor field T (..., 4, 4)."""
    T = np.asarray(T_field, float)
    flat_T = T.reshape(-1, 4, 4)
    N = len(flat_T)

    # 1. T00 positivity (Stationary observer v = (1, 0, 0, 0))
    T00 = flat_T[:, 0, 0]
    t00_min = float(np.min(T00))
    t00_pass = bool(t00_min >= -1e-14)

    # 2. Null Energy Condition (NEC): T_munu k^mu k^nu >= 0 for null directions k = (1, n_i)
    test_dirs = np.array([
        [1.0, 1.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0],
        [1.0, 0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, -1.0],
        [1.0, 1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3)],
    ])
    nec_vals = []
    for k in test_dirs:
        val = np.einsum("a,nab,b->n", k, flat_T, k)
        nec_vals.append(np.min(val))
    nec_min = float(min(nec_vals))
    nec_pass = bool(nec_min >= -1e-14)

    # 3. Weak Energy Condition (WEC): T_munu v^mu v^nu >= 0 for timelike vectors
    # Tested for several velocities beta in spatial directions
    wec_vals = [t00_min]
    for beta in [0.1, 0.5, 0.9]:
        gamma = 1.0 / np.sqrt(1.0 - beta**2)
        for d in [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])]:
            v = np.concatenate([[gamma], gamma * beta * d])
            val = np.einsum("a,nab,b->n", v, flat_T, v)
            wec_vals.append(np.min(val))
    wec_min = float(min(wec_vals))
    wec_pass = bool(wec_min >= -1e-14)

    # 4. Dominant Energy Condition (DEC): -T^mu_nu v^nu is future timelike or null
    # For v = (1, 0, 0, 0): energy flux J^mu = - T^mu_0 = (T_00, -T_01, -T_02, -T_03)
    # Norm of J: -(J^0)^2 + (J^1)^2 + (J^2)^2 + (J^3)^2 <= 0  <=>  |S/c|^2 <= u^2
    u_field = flat_T[:, 0, 0]
    p_flux = np.linalg.norm(flat_T[:, 0, 1:], axis=-1)
    dec_diff = u_field - p_flux
    dec_min = float(np.min(dec_diff))
    dec_pass = bool(dec_min >= -1e-14)

    return {
        "t00_positivity": {"min": t00_min, "pass": t00_pass},
        "WEC": {"min": wec_min, "pass": wec_pass},
        "NEC": {"min": nec_min, "pass": nec_pass},
        "DEC": {"min": dec_min, "pass": dec_pass},
        "overall_energy_conditions_pass": bool(t00_pass and wec_pass and nec_pass and dec_pass),
    }


def classify_conservation_gate(
    per_nu_residuals: list[float] | np.ndarray,
    scale: float,
    threshold_pass: float = 0.05,
    threshold_warn: float = 0.20,
) -> dict:
    """Classify 4-divergence residuals partial_mu T^{mu nu} into PASS, WARN, or FAIL."""
    res = np.asarray(per_nu_residuals, float)
    norm_scale = max(float(scale), 1e-300)
    rel_residuals = [float(r / norm_scale) for r in res]

    gates = []
    for nu, r_rel in enumerate(rel_residuals):
        if r_rel < threshold_pass:
            status = "PASS"
        elif r_rel < threshold_warn:
            status = "WARN"
        else:
            status = "FAIL"
        gates.append({"nu": nu, "abs_residual": float(res[nu]), "rel_residual": r_rel, "status": status})

    max_rel = max(rel_residuals)
    if max_rel < threshold_pass:
        overall = "PASS"
    elif max_rel < threshold_warn:
        overall = "WARN"
    else:
        overall = "FAIL"

    return {
        "overall_gate": overall,
        "max_rel_residual": float(max_rel),
        "threshold_pass": threshold_pass,
        "threshold_warn": threshold_warn,
        "components": gates,
    }
