"""Kerr + SC + condensate matter IC (v6.3-alpha2).

Combines:
- Kerr-Schild metric with smoothed M_eff(r) (kerr_sc_ic)
- SC condensate field (sc_condensate) normalized to perturbative amp
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .kerr_sc_ic import kerr_sc_state
from .sc_condensate import sc_condensate_profile, sc_condensate_potential


def kerr_sc_full_state(N, L, M=1.0, a=0.0, rc=1.0, rho_core=1e19,
                        m_pair_MeV=10.0, g_int=0.1, phi_amp_target=1e-4):
    """Return (state, phi_field, Pi_field, V_func)."""
    s = kerr_sc_state(N, L, M=M, a=a, rc=rc)
    phi_f, Pi_f, sigma, omega, A_amp = sc_condensate_profile(
        N=N, L=L, M=M, rho_core=rho_core, m_pair_MeV=m_pair_MeV)
    maxv = np.max(np.abs(phi_f))
    if maxv > 0:
        phi_f = phi_f * (phi_amp_target / maxv)
    V_func = sc_condensate_potential(m_pair_MeV=m_pair_MeV, g_int=g_int)
    return s, phi_f, Pi_f, V_func


__all__ = ["kerr_sc_full_state"]
