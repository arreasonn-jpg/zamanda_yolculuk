"""Z4c-style constraint damping for v2.0-alpha4.

Adds auxiliary field Theta and damping parameter kappa1, kappa2 to the
BSSN evolution. In the weak-damping limit, the correction to dK is
-kappa1 * alpha * H (Hamiltonian) and to dGamma is -kappa2 * M
(momentum). This is the minimal Z4c-inspired damping layer.
"""
from __future__ import annotations
import numpy as np


def damping_dK(alpha, H, kappa1=5.0):
    """Z4c dK correction: -kappa1 * alpha * H."""
    return -kappa1 * alpha * H


def damping_dGamma(M, kappa2=5.0):
    """Z4c dGamma correction: -kappa2 * M^i."""
    return -kappa2 * np.asarray(M, float)


def damping_energy(H, M):
    """Diagnostic: total constraint violation norm."""
    return float(abs(H) + np.linalg.norm(M))


__all__ = ["damping_dK", "damping_dGamma", "damping_energy"]
