"""Moving-puncture gauge (Bruegmann 2005):

  d_t alpha = -2 alpha K                          (1+log, moving)
  d_t beta^i = (3/4) B^i
  d_t B^i = d_t Gamma^i - eta B^i
"""
from __future__ import annotations
import numpy as np


def dt_alpha_mp(alpha, K, beta_dot_alpha=0.0):
    """d_t alpha = -2 alpha K + beta^i d_i alpha."""
    return -2.0 * alpha * K + beta_dot_alpha


def dt_beta_mp(B):
    """d_t beta^i = (3/4) B^i."""
    return 0.75 * np.asarray(B, float)


def dt_B_mp(dGamma, B, eta=1.0):
    """d_t B^i = dGamma^i - eta B^i."""
    return np.asarray(dGamma, float) - eta * np.asarray(B, float)


__all__ = ["dt_alpha_mp", "dt_beta_mp", "dt_B_mp"]
