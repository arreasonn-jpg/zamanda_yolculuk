"""Kerr + PTUM initial data (v6.2-alpha1).

Combines Kerr-Schild Cartesian metric with PTUM tachyon fields.

Kerr metric from bssn_ic_kerr_full.schwarzschild_puncture_state pattern
(via kerr_schild_state for KS coordinates).
PTUM fields: fuel profile Phi + equilibrium psi.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_ic_kerr_full import kerr_schild_state
from .ptum_field import fuel_profile, psi_equilibrium, m_eff_squared
from .matter_ptum import ptum_sources


def kerr_ptum_state(N, L, M=1.0, a=0.0, A_fuel=3.0, sigma=None,
                    m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Return (BSSNState, Phi, psi, Pi_Phi, Pi_psi, V_func)."""
    if abs(a) > M:
        raise ValueError("|a| <= M required")
    s = kerr_schild_state(N, L, M=M, a=a)
    if sigma is None:
        sigma = L / 6.0
    Phi = fuel_profile(N, L, Phi0=A_fuel, sigma=sigma)
    psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
    Z = np.zeros((N, N, N))
    return s, Phi, psi, Z.copy(), Z.copy()


def kerr_ptum_diagnostics(state, Phi, Pi_Phi, psi, Pi_psi, h,
                          V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0,
                          Eth=1.0, lam_psi=1.0):
    """Return dict with rho_min, n_neg, Kmax, and activation count."""
    src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, state.gtilde, state.phi, h,
                       V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                       lam_psi=lam_psi)
    meff = m_eff_squared(Phi, m2, alpha_p, Eth)
    return {
        "rho_min": float(src["rho"].min()),
        "rho_max": float(src["rho"].max()),
        "n_neg": int(np.sum(src["rho"] < 0)),
        "n_active": int(np.sum(meff < 0)),
        "Kmax": float(np.max(np.abs(state.K))),
        "alpha_min": float(state.alpha.min()),
    }


__all__ = ["kerr_ptum_state", "kerr_ptum_diagnostics"]
