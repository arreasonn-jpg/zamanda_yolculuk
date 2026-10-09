"""Cross-check: in the linear regime, BSSN perturbation h_ij should
match the retarded linearized-GR solution order-of-magnitude."""
from __future__ import annotations
import numpy as np
from .retarded_gr import C, retarded_hbar


def hbar_to_metric(hbar):
    """h_munu = hbar_munu - 1/2 eta_munu hbar (Lorenz gauge)."""
    hb = np.asarray(hbar, float)
    tr = np.trace(hb)  # eta-trace
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])  # (-,+,+,+)
    return hb - 0.5 * eta * tr


def linear_regime_ratio(h_bssn_amp, h_linear_amp):
    """Ratio nonlinear/linear amplitude. < 1 means safe linear regime."""
    return float(abs(h_bssn_amp) / max(abs(h_linear_amp), 1e-300))


__all__ = ["hbar_to_metric", "linear_regime_ratio"]
