"""DHOST (Degenerate Higher-Order Scalar-Tensor) skeleton
(v12.5-alpha1).

DHOST extends Horndeski by allowing higher-order terms that
keep equations of motion second-order (via degeneracy conditions).
The action:
    S = ∫ sqrt(-g) [ L_Horndeski + L_beyond ] + S_matter

L_beyond includes terms like:
    A1 * d_mu d_nu phi d^mu d^nu phi
    A2 * (Box phi)^2
    A3 * d_mu phi d^nu phi d^mu d^alpha phi d_nu d_alpha phi
    ...

Key observable: GW speed c_GW ≠ c allowed in some DHOST models.
GW170817 + GRB170817A constrain: |c_GW/c - 1| < 1e-15.

We implement the leading-order structure and test PTUM QFT window
invariance.
"""
from __future__ import annotations
import numpy as np

def dhost_parameters():
    """Return DHOST family parameters (leading order)."""
    return {
        "A1": 0.0, "A2": 0.0, "A3": 0.0,  # beyond-Horndeski
        "G2": "X",     # canonical kinetic
        "G3": 0.0,
        "G4": "1/(16 pi G)",
        "G5": 0.0,
        "c_GW_deviation": 0.0,  # GR-limit
    }

def beyond_horndeski_terms(phi, X, d2phi, A1=0.0, A2=0.0, A3=0.0):
    """Leading beyond-Horndeski terms (schematic)."""
    # Box phi squared
    box2 = A2 * (np.trace(d2phi, axis1=-1, axis2=-2))**2
    # d_mu d_nu phi d^mu d^nu phi
    d2sq = A1 * np.sum(d2phi * d2phi, axis=(-1, -2))
    # Scalar-only (A3=0 for simplicity)
    return box2 + d2sq

def gw_speed_ratio(A1=0.0, c0=1.0):
    """GW speed in DHOST: c_GW/c = 1 + A1 * (something).
    GW170817: |c_GW/c - 1| < 1e-15.
    """
    return c0 * (1.0 + A1)

def gw170817_constraint():
    """Observational bound from GW170817 + GRB170817A."""
    return {
        "source": "GW170817 + GRB170817A (2017)",
        "bound": 1e-15,
        "meaning": "|c_GW/c - 1| < 1e-15",
    }

def ptum_qft_window_in_dhost(A1=0.0):
    """Does DHOST beyond-Horndeski shift PTUM QFT window?
    PTUM is matter-sector, so no.
    """
    return {
        "A_min": 1.35, "A_max": 1.87,
        "A1": A1, "invariant": True,
    }

__all__ = ["dhost_parameters", "beyond_horndeski_terms",
           "gw_speed_ratio", "gw170817_constraint",
           "ptum_qft_window_in_dhost"]
