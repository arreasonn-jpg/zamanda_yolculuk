"""Quantum inequality (Ford-Roman) diagnostic (v6.0-alpha8).

Ford-Roman bound for massless scalar field (4D):
  ∫ rho(t) g(t) dt  >=  -3 / (32 pi^2 tau^4)
for a sampling function g(t) = tau / (pi (t^2 + tau^2)) with
∫ g dt = 1.

Test: compute the time-integrated negative energy density sampled
by g(t) over a sequence of snapshots, and compare to the bound.

This is a diagnostic, NOT a proof. A PTUM discharge that violates
the Ford-Roman bound would indicate a genuine conflict with quantum
field theory.
"""
from __future__ import annotations
import numpy as np

HBAR = 1.054571817e-34
C = 2.99792458e8


def sampling_function(t, tau):
    """Lorentzian sampling: g(t) = tau / (pi (t^2 + tau^2))."""
    return tau / (np.pi * (t**2 + tau**2))


def ford_roman_bound(tau):
    """Massless scalar 4D: -3 / (32 pi^2 tau^4)."""
    return -3.0 / (32.0 * np.pi**2 * tau**4)


def quantum_inequality_value(rho_history, t_history, tau):
    """Compute integral of rho(t) g(t) dt over the sample."""
    rho = np.asarray(rho_history, float)
    t = np.asarray(t_history, float)
    g = sampling_function(t, tau)
    # Trapezoid integration
    integrand = rho * g
    val = np.trapezoid(integrand, t)
    return float(val)


def check_ford_roman(rho_history, t_history, tau,
                     dimension=4, N_field_components=1):
    """Return dict with inequality value and bound.

    N_field_components: number of scalar fields contributing
    (multiplies bound by N).
    """
    val = quantum_inequality_value(rho_history, t_history, tau)
    bound = N_field_components * ford_roman_bound(tau)
    return {
        "value": val,
        "bound": bound,
        "violates": val < bound,
        "ratio": val / bound if bound != 0 else np.inf,
        "tau": tau,
    }


__all__ = ["sampling_function", "ford_roman_bound",
           "quantum_inequality_value", "check_ford_roman"]
