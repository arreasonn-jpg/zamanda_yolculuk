"""Starobinsky f(R) = R + alpha R^2 as scalar-tensor (v8.2-alpha1).

Equivalent scalar-tensor form (in Planck units):
  phi = sqrt(3/2) ln(1 + 2 alpha R)
  V(phi) = (1/(8 alpha)) (1 - e^{-sqrt(2/3) phi})^2

This is the R^2 inflaton potential.
"""
from __future__ import annotations
import numpy as np


def V_starobinsky(phi, alpha=1.0):
    """Starobinsky potential V(phi)."""
    phi = np.asarray(phi, float)
    pref = 1.0 / (8.0 * alpha)
    return pref * (1.0 - np.exp(-np.sqrt(2.0/3.0) * phi))**2


def dV_starobinsky(phi, alpha=1.0):
    phi = np.asarray(phi, float)
    pref = 1.0 / (8.0 * alpha)
    s = np.sqrt(2.0/3.0)
    e = np.exp(-s * phi)
    return pref * 2.0 * (1.0 - e) * s * e


def R_from_phi(phi, alpha=1.0):
    """Inverse: R(phi) = (e^{sqrt(2/3) phi} - 1) / (2 alpha)."""
    s = np.sqrt(2.0/3.0)
    return (np.exp(s * np.asarray(phi, float)) - 1.0) / (2.0 * alpha)


def slow_roll_observables(alpha, N_efolds=55):
    """Approximate n_s and r for Starobinsky."""
    # Standard results at N=55: n_s ~ 0.965, r ~ 0.004 for alpha=1
    ns = 1.0 - 2.0 / N_efolds - 3.0 / (2.0 * N_efolds**2)
    r = 12.0 / N_efolds**2
    return {"n_s": ns, "r": r, "N_efolds": N_efolds}


__all__ = ["V_starobinsky", "dV_starobinsky", "R_from_phi",
           "slow_roll_observables"]
