"""Teleparallel Equivalent of GR (TEGR) skeleton (v12.5-alpha2).

Instead of curvature (R), teleparallel gravity uses torsion (T).
Equivalent to GR at the level of equations of motion, but the
underlying geometry is different (Weitzenbock connection).

Action: S = (1/(16 pi G)) ∫ e T d^4x + S_matter
where e = det(e^a_mu) and T is the torsion scalar.

Modified: f(T) gravity is the teleparallel analog of f(R).
Observable differences arise at higher order.
"""
from __future__ import annotations
import numpy as np

def tetrad_flat():
    """Flat-space tetrad e^a_mu = delta^a_mu."""
    return np.eye(4)

def torsion_scalar_flat():
    """Torsion scalar for flat tetrad."""
    return 0.0

def f_T_gravity(T, alpha=0.0):
    """f(T) = T + alpha T^2 (Starobinsky-like in teleparallel)."""
    return T + alpha * T**2

def f_T_GR_limit():
    """alpha=0 -> TEGR -> GR."""
    return {"alpha": 0.0, "equivalent_to_GR": True}

def ptum_qft_window_in_fT(alpha=0.0):
    """PTUM QFT window: matter-sector, invariant."""
    return {
        "A_min": 1.35, "A_max": 1.87,
        "alpha": alpha, "invariant": True,
    }

__all__ = ["tetrad_flat", "torsion_scalar_flat",
           "f_T_gravity", "f_T_GR_limit", "ptum_qft_window_in_fT"]
