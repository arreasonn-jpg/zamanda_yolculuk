"""Schwarzschild Condensation regular core (v3.0-alpha4).

Reference: Esa (2026), "Schwarzschild Condensation" — QCD-inspired
BEC model for non-singular BH interiors. This module implements the
effective metric in geometric units (G = c = 1):

  Interior (r < R0): de Sitter-like
    f_in(r)  = 1 - (r/ell)^2,   ell^2 = 3 / (8 pi rho_0)
  Exterior (r > R0): Schwarzschild
    f_out(r) = 1 - 2 M / r

with matching at R0 = ell (where f_in = 0).

Curvature diagnostics:
  Ricci scalar R(r), Kretschmann K(r) for the static spherically
  symmetric metric ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2.
"""
from __future__ import annotations
import numpy as np


def schwarzschild_radius(M):
    """Rs = 2M in geometric units."""
    return 2.0 * M


def de_sitter_length(rho0):
    """ell^2 = 3 / (8 pi rho0) in geometric units."""
    return np.sqrt(3.0 / (8.0 * np.pi * rho0))


def f_SC(r, M, R0, ell):
    """Effective metric function f(r) = -g_tt = 1/g_rr."""
    r = np.asarray(r, float)
    Rs = schwarzschild_radius(M)
    f_in = 1.0 - (r / ell) ** 2
    f_out = 1.0 - Rs / np.maximum(r, 1e-30)
    return np.where(r < R0, f_in, f_out)


def rho_SC(r, M, R0, ell):
    """Effective energy density: constant rho0 inside R0, 0 outside.

    rho0 = 3 / (8 pi ell^2) from de Sitter relation.
    """
    rho0 = 3.0 / (8.0 * np.pi * ell**2)
    return np.where(np.asarray(r, float) < R0, rho0, 0.0)


def ricci_scalar_SC(r, M, R0, ell, h=1e-6):
    """Ricci scalar of the SC metric via finite differences.

    For static spherical metric:
      R = -f'' - (4/r) f' + 2(1-f)/r^2
    """
    r = np.asarray(r, float)
    f = f_SC(r, M, R0, ell)
    fp  = (f_SC(r + h, M, R0, ell) - f_SC(r - h, M, R0, ell)) / (2.0 * h)
    fpp = (f_SC(r + h, M, R0, ell) - 2.0 * f + f_SC(r - h, M, R0, ell)) / h**2
    return -fpp - 4.0 * fp / np.maximum(r, 1e-30) + 2.0 * (1.0 - f) / np.maximum(r, 1e-30)**2


def kretschmann_SC(r, M, R0, ell, h=1e-6):
    """Kretschmann-like proxy:
      K = (f'')^2 + 4 (f'/r)^2 + 4 (1-f)/r^2)^2
    """
    r = np.asarray(r, float)
    f = f_SC(r, M, R0, ell)
    fp  = (f_SC(r + h, M, R0, ell) - f_SC(r - h, M, R0, ell)) / (2.0 * h)
    fpp = (f_SC(r + h, M, R0, ell) - 2.0 * f + f_SC(r - h, M, R0, ell)) / h**2
    return fpp**2 + 4.0 * (fp / np.maximum(r, 1e-30))**2 + 4.0 * ((1.0 - f) / np.maximum(r, 1e-30)**2)**2


def junction_surface_density(M, R0, ell):
    """Surface density from Israel junction (thin-shell approximation).

    sigma_sh = -(1/(4 pi G)) ([K] - h [K_ab]) evaluated on shell.
    Simplified form for dS-Sch matching:
      sigma_sh = -(1/(2 pi R0)) (sqrt(f_out) - sqrt(f_in))
    """
    Rs = schwarzschild_radius(M)
    f_in_val = 1.0 - (R0 / ell) ** 2
    f_out_val = 1.0 - Rs / R0
    if f_in_val < 0 or f_out_val < 0:
        return np.nan
    return -(1.0 / (2.0 * np.pi * R0)) * (np.sqrt(f_out_val) - np.sqrt(f_in_val))


def core_parameters(M, rho0):
    """Return (Rs, R0, ell, delta) given M and interior density rho0 in geometric units."""
    Rs = schwarzschild_radius(M)
    ell = de_sitter_length(rho0)
    # Matching at R0 = ell where f_in(R0) = 0
    R0 = ell
    delta = R0 / Rs - 1.0
    return {"Rs": Rs, "R0": R0, "ell": ell, "delta": delta}


__all__ = [
    "schwarzschild_radius", "de_sitter_length", "f_SC", "rho_SC",
    "ricci_scalar_SC", "kretschmann_SC", "junction_surface_density",
    "core_parameters",
]
