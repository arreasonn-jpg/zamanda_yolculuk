"""QFT-OK window invariance under LQC (v8.2-alpha3).

Combines LQC Friedmann correction with PTUM discharge.
LQC correction to H^2 is (1 - rho/rho_c).
For PTUM rho ~ 1e19 kg/m^3 and rho_c ~ 1e96 kg/m^3,
the correction factor is 1 - 1e-77, i.e., zero to machine precision.
"""
from __future__ import annotations
import numpy as np
from .lqc_scalar_evolution import RHO_C

# Representative PTUM densities (from v7.4 four-regime analysis)
A_TO_RHO_MIN = {
    1.4: -0.137,
    1.6: -0.488,
    1.8: -1.105,
    2.0: -2.069,
}

# PTUM energy density in physical units (kg/m^3)
# PTUM code units: rho ~ 1 corresponds roughly to nuclear density
RHO_PTUM_UNIT = 1e19


def lqc_correction_factor(rho_physical, rho_c=RHO_C):
    """1 - rho/rho_c: how much LQC suppresses H^2."""
    return 1.0 - rho_physical / rho_c


def ptum_rho_physical(A):
    """PTUM energy density in physical units (kg/m^3)."""
    rho_code = A_TO_RHO_MIN.get(A)
    if rho_code is None:
        raise ValueError(f"A={A} not in table")
    return abs(rho_code) * RHO_PTUM_UNIT


def qft_window_shift_due_to_lqc(A, rho_c=RHO_C):
    """How much does LQC shift the PTUM QFT ratio?"""
    rho_phys = ptum_rho_physical(A)
    corr = lqc_correction_factor(rho_phys, rho_c)
    # QFT ratio scales with rho_min (linear law from v7.3)
    # LQC correction reduces rho by factor corr -> ratio scales by corr
    return {
        "A": A,
        "rho_ptum_kg_per_m3": rho_phys,
        "lqc_correction": corr,
        "shift": 1.0 - corr,  # relative shift
    }


__all__ = [
    "lqc_correction_factor", "ptum_rho_physical",
    "qft_window_shift_due_to_lqc", "A_TO_RHO_MIN", "RHO_PTUM_UNIT",
]
