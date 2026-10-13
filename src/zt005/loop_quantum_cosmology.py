"""Loop Quantum Cosmology (LQC) effective dynamics (v8.1-alpha1).

LQC replaces the classical Friedmann equation with a holonomy-corrected
version:
    H^2 = (8 pi G / 3) rho * (1 - rho / rho_c)

where rho_c ~ 0.41 rho_Planck is the critical density. At rho = rho_c,
H = 0 (bounce). For rho << rho_c, classical Friedmann is recovered.

We implement:
- LQC-modified Friedmann equation
- Bounce condition (H = 0 at rho_c)
- Effective scalar field evolution in LQC background
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34

# Planck density (kg/m^3): rho_P = c^5 / (hbar G^2)
RHO_PLANCK = C**5 / (HBAR * G**2)

# LQC critical density: rho_c ~ 0.41 rho_Planck (Ashtekar 2006)
RHO_C = 0.41 * RHO_PLANCK


def friedmann_lqc(rho, rho_c=RHO_C):
    """LQC-modified Friedmann: H^2 = (8 pi G / 3) rho (1 - rho / rho_c)."""
    return (8.0 * np.pi * G / 3.0) * rho * (1.0 - rho / rho_c)


def friedmann_classical(rho):
    """Classical Friedmann: H^2 = (8 pi G / 3) rho."""
    return (8.0 * np.pi * G / 3.0) * rho


def bounce_density(rho_c=RHO_C):
    """Density at which H = 0 (bounce)."""
    return rho_c


def is_bouncing(rho, rho_c=RHO_C):
    """True if rho > rho_c/2 (i.e., past the classical maximum)."""
    return rho > 0.5 * rho_c


def hubble_ratio(rho, rho_c=RHO_C):
    """H_LQC^2 / H_classical^2 = 1 - rho/rho_c."""
    if rho <= 0:
        return 1.0
    return 1.0 - rho / rho_c


__all__ = [
    "RHO_PLANCK", "RHO_C", "friedmann_lqc", "friedmann_classical",
    "bounce_density", "is_bouncing", "hubble_ratio",
]
