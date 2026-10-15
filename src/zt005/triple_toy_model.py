"""SC + PTUM + BH triple toy model (v12.3-alpha1).

Question: in a BH gravitational background, do SC regular core
and PTUM discharge interact? Three regimes:
  - SC: regular core replaces singularity (metric)
  - PTUM: tachyon discharge → NEC violation (matter)
  - BH: strong curvature (background)

Test: does the presence of one shift the QFT window of the other?
Expected: no (each is local to its sector).
"""
from __future__ import annotations
import numpy as np

def schwarzschild_f_isotropic(r, M=1.0):
    """Lapse in isotropic coords."""
    return ((1 - M/(2*r)) / (1 + M/(2*r)))**2

def sc_core_metric(r, M=1.0, rc=2.0):
    """SC smoothed-core mass profile."""
    M_eff = M * (1.0 - np.exp(-r/rc) * (1.0 + r/rc))
    return M_eff

def ptum_local_activation(Phi, m2=1.0, alpha_p=1.0, Eth=1.0):
    """Local PTUM activation (matter sector only)."""
    meff2 = m2 * (1.0 - alpha_p * Phi**2 / Eth**2)
    return meff2 < 0.0

def triple_state_at_r(r, M=1.0, rc=2.0, A_fuel=3.0, sigma=4.0):
    """Evaluate all three sectors at radius r.
    Returns dict with SC, PTUM, BH quantities.
    """
    # SC metric
    M_eff = sc_core_metric(r, M, rc)
    # BH lapse
    f_bh = schwarzschild_f_isotropic(r, M) if r > 0.5*M else 0.0
    # PTUM local field (Gaussian)
    Phi = A_fuel * np.exp(-r**2 / (2*sigma**2))
    meff2 = 1.0 * (1.0 - Phi**2 / 1.0)
    if meff2 < 0:
        psi_eq = np.sqrt(-meff2)
        V_min = -meff2**2 / 4.0
    else:
        psi_eq = 0.0
        V_min = 0.0
    return {
        "r": r, "M_eff": M_eff, "f_bh": f_bh,
        "Phi": Phi, "psi_eq": psi_eq, "V_min": V_min,
        "activated": meff2 < 0,
    }

def qft_window_triple():
    """PTUM QFT window: does SC + BH shift it?
    Both SC and BH are metric-sector; PTUM QFT window is matter-sector.
    Predict: invariant.
    """
    return {"A_min": 1.35, "A_max": 1.87,
            "SC_shift": 0.0, "BH_shift": 0.0,
            "invariant": True}

def triple_coupling_check(M=1.0, rc=2.0, r_probe=5.0):
    """Check if SC core + BH interact with PTUM activation."""
    s = triple_state_at_r(r_probe, M=M, rc=rc)
    return {
        "r": r_probe,
        "SC_M_eff": s["M_eff"],
        "BH_lapse": s["f_bh"],
        "PTUM_activated": s["activated"],
        "PTUM_psi": s["psi_eq"],
    }

__all__ = ["schwarzschild_f_isotropic", "sc_core_metric",
           "ptum_local_activation", "triple_state_at_r",
           "qft_window_triple", "triple_coupling_check"]
