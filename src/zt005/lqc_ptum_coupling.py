"""LQC + PTUM coupling (v8.2-alpha2).

LQC modifies the Friedmann equation; PTUM provides the matter.
We test whether LQC's rho_c bound affects the PTUM QFT-OK window.

Physics: LQC is a global (background) effect; PTUM is local.
Expectation: PTUM QFT threshold unchanged (as with Kerr/BD).
"""
from __future__ import annotations
import numpy as np
from .lqc_scalar_evolution import H_squared_lqc, RHO_C
from .matter.ptum_field import V_fuel, m_eff_squared


def combined_density_lqc_ptum(Phi, Pi_Phi, psi, Pi_psi,
                                V0=1.0, lam=0.05, m2=1.0,
                                alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Total energy density (homogeneous) for LQC + PTUM."""
    V_f = V_fuel(Phi, V0=V0, lam=lam)
    meff2 = m_eff_squared(Phi, m2=m2, alpha=alpha_p, Eth=Eth)
    V_p = 0.5 * meff2 * psi**2 + 0.25 * lam_psi * psi**4
    rho = 0.5*Pi_Phi**2 + 0.5*Pi_psi**2 + V_f + V_p
    return float(rho)


def is_below_rho_c(rho, rho_c=RHO_C):
    """True if rho << rho_c (classical regime)."""
    return rho < 0.01 * rho_c


__all__ = ["combined_density_lqc_ptum", "is_below_rho_c"]
