"""Nonlinear GR skeleton: scalar source, second-order back-reaction,
linear->nonlinear crossover. Minimal; not a BSSN solver."""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

C = 2.99792458e8
G = 6.67430e-11


def scalar_field_stress_energy(phi, dphi, g, ginv, V=0.0):
    """T_munu = d_mu phi d_nu phi - g_munu (1/2 g^ab d_a phi d_b phi + V)."""
    phi = np.asarray(phi)
    dphi = np.asarray(dphi, dtype=float)
    g = np.asarray(g, dtype=float)
    ginv = np.asarray(ginv, dtype=float)
    if dphi.ndim == 1:
        K = 0.5 * float(dphi @ ginv @ dphi)
        Vv = float(V(phi)) if callable(V) else float(V)
        return np.outer(dphi, dphi) - g * (K + Vv)
    K = 0.5 * np.einsum('ni,nij,nj->n', dphi, ginv, dphi)
    Vv = np.array([V(p) for p in phi]) if callable(V) else np.full(dphi.shape[0], float(V))
    return np.einsum('ni,nj->nij', dphi, dphi) - g * (K + Vv)[:, None, None]


def nonlinearity_parameter(h):
    """epsilon = |h_munu|_max; nonlinear regime when epsilon ~ O(1)."""
    return float(np.max(np.abs(h)))


def orders_to_nonlinear(h_current, h_nonlinear=1.0):
    """log10(h_nl / h_current) — orders of magnitude to close."""
    return float(np.log10(h_nonlinear / max(abs(h_current), 1e-300)))


def second_order_backreaction(h_amp, k_wave):
    """Landau-Lifshitz pseudotensor order-of-magnitude estimate.
    Linear: (c^4/G) k^2 h^2 ; Nonlinear (2nd order): (c^4/G) k^2 h^4."""
    lin = (C**4 / G) * k_wave**2 * h_amp**2
    nl  = (C**4 / G) * k_wave**2 * h_amp**4
    return {"linear_J_per_m3": lin, "nonlinear_J_per_m3": nl,
            "ratio": nl / max(lin, 1e-300)}


def crossover_h():
    """h at which linear and nonlinear terms become comparable."""
    return 1.0


@dataclass
class NonlinearBudget:
    h_current: float = 1e-46
    h_nonlinear: float = 1.0
    def summary(self) -> dict:
        return {
            "h_current": self.h_current,
            "h_nonlinear": self.h_nonlinear,
            "epsilon_current": self.h_current,
            "orders_to_nonlinear": orders_to_nonlinear(self.h_current, self.h_nonlinear),
        }


__all__ = ["scalar_field_stress_energy", "nonlinearity_parameter",
           "orders_to_nonlinear", "second_order_backreaction",
           "crossover_h", "NonlinearBudget"]
