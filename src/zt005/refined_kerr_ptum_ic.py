"""Refined Kerr + PTUM IC helper (v6.2-alpha2).

Builds coarse + fine Kerr-Schild PTUM states:
- Coarse: wide domain L_c, N_c
- Fine: narrow domain L_f centered, N_f
Both use same physical M, a, A_fuel, sigma (in physical units).
"""
from __future__ import annotations
import numpy as np
from .kerr_ptum_ic import kerr_ptum_state


def refined_kerr_ptum_ic(N_c, L_c, N_f, L_f, M=1.0, a=0.0, A_fuel=3.0,
                         sigma_phys=4.0, m2=1.0, alpha_p=1.0,
                         Eth=1.0, lam_psi=1.0):
    """Return (s_c, Phi_c, psi_c, s_f, Phi_f, psi_f)."""
    s_c, Phi_c, psi_c, P1c, P2c = kerr_ptum_state(
        N=N_c, L=L_c, M=M, a=a, A_fuel=A_fuel, sigma=sigma_phys,
        m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    s_f, Phi_f, psi_f, P1f, P2f = kerr_ptum_state(
        N=N_f, L=L_f, M=M, a=a, A_fuel=A_fuel, sigma=sigma_phys,
        m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    return s_c, Phi_c, psi_c, s_f, Phi_f, psi_f


__all__ = ["refined_kerr_ptum_ic"]
