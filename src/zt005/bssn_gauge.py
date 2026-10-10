"""Gauge evolution for BSSN v2.0-alpha3.

1+log lapse:        d_t alpha = -2 alpha K          (+ L_beta alpha)
Gamma-driver shift: d_t beta^i = (3/4) B^i
                    d_t B^i    = d_t Gamma^i - eta B^i
B^i stored as auxiliary; we use a simplified form with B^i folded.
"""
from __future__ import annotations
import numpy as np


def dt_alpha_1pluslog(alpha, K, beta_dot_alpha=0.0):
    """d_t alpha = -2 alpha K + beta^i d_i alpha  (shift term passed in)."""
    return -2.0 * alpha * K + beta_dot_alpha


def dt_beta_gamma_driver(dt_Gamma, beta, eta=1.0):
    """d_t beta^i = (3/4) (dt_Gamma^i - eta beta^i)."""
    return 0.75 * (np.asarray(dt_Gamma, float) - eta * np.asarray(beta, float))


__all__ = ["dt_alpha_1pluslog", "dt_beta_gamma_driver"]
