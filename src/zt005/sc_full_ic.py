"""Full SC integration: SC metric + SC condensate matter (v4.0-alpha3).

Combines:
- SC smoothed-core metric (sc_like_ic): psi(r) = 1 + M_eff(r)/(2r)
- SC condensate scalar field (sc_condensate): Gaussian, normalized

The scalar field is rescaled so that its stress-energy stays perturbative
relative to the SC metric (ratio ~ 1e-3).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .sc_like_ic import sc_like_state
from .matter.sc_condensate import sc_condensate_profile, sc_condensate_potential


def sc_full_state(N, L, M=1.0, rc=1.0, rho_core=1e19, m_pair_MeV=10.0,
                  g_int=0.1, phi_amp_target=1e-4):
    """SC metric + SC condensate scalar field.

    The scalar amplitude is rescaled to phi_amp_target (default 1e-4),
    small enough that its T_mu_nu is a perturbation of the metric.
    """
    metric_state = sc_like_state(N, L, M=M, rc=rc)
    phi_f, Pi_f, sigma, omega, A_amp = sc_condensate_profile(
        N=N, L=L, M=M, rho_core=rho_core, m_pair_MeV=m_pair_MeV)
    # normalize to target amplitude
    maxv = np.max(np.abs(phi_f))
    if maxv > 0:
        phi_f = phi_f * (phi_amp_target / maxv)
    V_func = sc_condensate_potential(m_pair_MeV=m_pair_MeV, g_int=g_int)
    return metric_state, phi_f, Pi_f, V_func


__all__ = ["sc_full_state"]
