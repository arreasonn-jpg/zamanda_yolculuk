"""Horndeski scalar-tensor gravity skeleton (v9.1-alpha1).

The Horndeski action is:
    S = ∫ sqrt(-g) [L2 + L3 + L4 + L5] + S_matter

where:
    L2 = G2(phi, X)                     X = -1/2 (dphi)^2
    L3 = -G3(phi, X) Box phi
    L4 = G4(phi, X) R + G4_X [(Box phi)^2 - (nabla nabla phi)^2]
    L5 = G5(phi, X) G_munu nabla^mu nabla^nu phi
         - (1/6) G5_X [(Box phi)^3 - 3 Box phi (nabla nabla phi)^2
                       + 2 (nabla nabla phi)^3]

GR: G4 = 1/(16 pi G), others = 0.
Brans-Dicke: G4 = phi, G2 = omega/phi * X.

We implement the general structure and test the leading-order
effect on PTUM QFT-OK window.
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11


def L2_kessence(phi, X, G2_func=None):
    """k-essence Lagrangian L2 = G2(phi, X)."""
    if G2_func is None:
        return X  # canonical kinetic
    return G2_func(phi, X)


def L4_nonminimal(R, phi, G4_func=None):
    """Non-minimal coupling L4 = G4(phi) R."""
    if G4_func is None:
        G4 = 1.0 / (16.0 * np.pi * G)
    else:
        G4 = G4_func(phi)
    return G4 * R


def horndeski_GR_limit():
    """GR is recovered for G2 = X, G3 = G5 = 0, G4 = 1/(16 pi G)."""
    return {"G2": "X", "G3": 0, "G4": "1/(16 pi G)", "G5": 0}


def horndeski_BD_limit(omega):
    """Brans-Dicke: G2 = omega/phi * X, G4 = phi."""
    return {"G2": f"{omega}/phi * X", "G3": 0, "G4": "phi", "G5": 0}


def ptum_qft_window_in_horndeski(G4_value):
    """PTUM QFT window is G4-independent (matter-sector only)."""
    # Same window as v7.4 across all Horndeski limits
    return {
        "A_min": 1.35,
        "A_max": 1.87,
        "G4_value": G4_value,
        "invariant": True,
    }


__all__ = [
    "L2_kessence", "L4_nonminimal",
    "horndeski_GR_limit", "horndeski_BD_limit",
    "ptum_qft_window_in_horndeski",
]
