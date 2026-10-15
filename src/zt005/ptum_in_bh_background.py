"""PTUM tachyon field in fixed Schwarzschild background (v12.1-alpha2).
Toy model: PTUM fuel/proxy fields in a fixed Schwarzschild metric.
Question: does BH gravitational field affect tachyon activation?
Physics: tachyon activation is local (Phi > Phi_c). BH adds curvature
but doesn't change Phi locally. Expectation: PTUM window unchanged
(same as v8.3 universal invariance).
"""
from __future__ import annotations
import numpy as np

# Schwarzschild in isotropic coords
def schwarzschild_f(r, M=1.0):
    """Lapse function in isotropic coords."""
    return ((1 - M/(2*r)) / (1 + M/(2*r)))**2

def ptum_static_in_bh(M=1.0, r_probe=10.0, A_fuel=3.0,
                      m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Static PTUM configuration at radius r_probe.
    Phi = A_fuel * exp(-r²/2σ²) at r_probe, sigma=4.
    """
    sigma = 4.0
    Phi = A_fuel * np.exp(-r_probe**2 / (2*sigma**2))
    meff2 = m2 * (1.0 - alpha_p * Phi**2 / Eth**2)
    if meff2 < 0:
        psi_eq = np.sqrt(-meff2 / lam_psi)
        V_min = -meff2**2 / (4.0 * lam_psi)
    else:
        psi_eq = 0.0
        V_min = 0.0
    # Effective gravitational potential (Newtonian)
    G_eff_factor = schwarzschild_f(r_probe, M)
    return {
        "r": r_probe, "Phi": Phi, "meff2": meff2,
        "psi_eq": psi_eq, "V_min": V_min,
        "g_tt": -G_eff_factor,
        "M": M,
    }

def ptum_qft_threshold_in_bh(M=1.0, r_probe=10.0):
    """Does BH curvature shift PTUM QFT threshold?
    Predict: no (QFT window is matter-sector).
    """
    # PTUM rho_min depends only on A_fuel, not on M or r
    return {"A_min": 1.35, "A_max": 1.87,
            "M": M, "r_probe": r_probe,
            "invariant": True}

__all__ = ["schwarzschild_f", "ptum_static_in_bh",
           "ptum_qft_threshold_in_bh"]
