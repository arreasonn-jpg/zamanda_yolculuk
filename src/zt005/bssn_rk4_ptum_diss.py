"""RK4 PTUM driver with Kreiss-Oliger dissipation on PTUM fields
(v5.0-alpha4). Adds eps_KO dissipation to Phi, Pi_Phi, psi, Pi_psi."""
from __future__ import annotations
import numpy as np
from .bssn_rk4_ptum import rk4_step_ptum
from .bssn_dissipation import ko_dissipation


def rk4_step_ptum_diss(state, Phi, Pi_Phi, psi, Pi_psi, h, dt,
                       V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0,
                       lam_psi=1.0, eps_KO=0.1):
    """RK4 step with KO dissipation applied to PTUM fields after each step."""
    s2, Phi2, Pi_Phi2, psi2, Pi_psi2 = rk4_step_ptum(
        state, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    scale = dt / h
    Phi2 = Phi2 + scale * eps_KO * ko_dissipation(Phi2, h)
    Pi_Phi2 = Pi_Phi2 + scale * eps_KO * ko_dissipation(Pi_Phi2, h)
    psi2 = psi2 + scale * eps_KO * ko_dissipation(psi2, h)
    Pi_psi2 = Pi_psi2 + scale * eps_KO * ko_dissipation(Pi_psi2, h)
    return s2, Phi2, Pi_Phi2, psi2, Pi_psi2


__all__ = ["rk4_step_ptum_diss"]
