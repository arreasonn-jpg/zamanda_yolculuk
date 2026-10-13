"""Hu-Sawicki f(R) dark energy model (v8.3-alpha1).

f(R) = R - m^2 c1 (R/m^2)^n / (c2 (R/m^2)^n + 1)

- High curvature: f(R) -> R - 2 Lambda_eff (GR + cosmological constant)
- Low curvature: modified gravity regime
- Chameleon screening: solar system constraints satisfied

Parameters typically n=1, c1/c2^2 ~ Lambda/m^2.
"""
from __future__ import annotations
import numpy as np

# Cosmological parameters (Planck 2018)
H0 = 2.18e-18  # s^-1
OMEGA_M = 0.315
OMEGA_L = 0.685
RHO_C_CRIT = 8.6e-27  # kg/m^3 (critical density today)
LAMBDA_EFF = 1.1e-52  # m^-2 (cosmological constant)


def f_R(R, n=1.0, c1=None, c2=None, m2=None):
    """Hu-Sawicki f(R)."""
    if m2 is None:
        m2 = LAMBDA_EFF
    if c1 is None:
        c1 = 1.0
    if c2 is None:
        c2 = 1.0
    x = R / m2
    return R - m2 * c1 * x**n / (c2 * x**n + 1.0)


def f_R_deriv(R, n=1.0, c1=None, c2=None, m2=None):
    """df/dR for Hu-Sawicki."""
    if m2 is None:
        m2 = LAMBDA_EFF
    if c1 is None:
        c1 = 1.0
    if c2 is None:
        c2 = 1.0
    x = R / m2
    xp = c2 * x**n + 1.0
    # f = R - m2 * c1 * x^n / xp
    # df/dR = 1 - c1 * n * x^{n-1} / xp + c1 * c2 * n * x^{2n-1} / xp^2
    return 1.0 - c1 * n * x**(n-1) / xp + c1 * c2 * n * x**(2*n-1) / xp**2


def screening_ratio(R_bg, R_screen):
    """f_R(R_bg) / f_R(R_screen) - measures chameleon suppression."""
    return abs(f_R_deriv(R_bg)) / max(abs(f_R_deriv(R_screen)), 1e-300)


__all__ = ["f_R", "f_R_deriv", "screening_ratio",
           "H0", "OMEGA_M", "OMEGA_L", "LAMBDA_EFF"]
