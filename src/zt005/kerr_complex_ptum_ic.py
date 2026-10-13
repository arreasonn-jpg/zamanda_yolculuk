"""Kerr metric + complex scalar + PTUM IC (v7.1-alpha2).

Combines:
- Kerr-Schild metric (kerr_schild_state from bssn_ic_kerr_full)
- Complex scalar Psi = Psi_R + i Psi_I (perturbative)
- PTUM fields (Phi, psi) at equilibrium
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_ic_kerr_full import kerr_schild_state
from .ptum_field import fuel_profile, psi_equilibrium_smooth


def kerr_complex_ptum_state(N, L, M=1.0, a=0.0,
                              # complex scalar
                              Psi_R_amp=1e-4, Psi_I_amp=0.5e-4,
                              sigma_Psi=None,
                              # PTUM
                              A_fuel=3.0, sigma_fuel=None,
                              m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Return dict with BSSN state + all matter fields."""
    if abs(a) > M:
        raise ValueError("|a| <= M required")
    s = kerr_schild_state(N, L, M=M, a=a)
    if sigma_Psi is None:
        sigma_Psi = L / 10.0
    if sigma_fuel is None:
        sigma_fuel = L / 5.0

    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2

    # Complex scalar (perturbative Gaussian)
    PR = Psi_R_amp * np.exp(-R2 / (2.0 * sigma_Psi**2))
    PI = Psi_I_amp * np.exp(-R2 / (2.0 * sigma_Psi**2))
    Zc = np.zeros((N, N, N))

    # PTUM
    Phi = fuel_profile(N, L, Phi0=A_fuel, sigma=sigma_fuel)
    psi = psi_equilibrium_smooth(Phi, m2=m2, alpha=alpha_p, Eth=Eth,
                                   lam_psi=lam_psi)

    return {
        "state": s,
        "PR": PR, "PI": PI, "PiR": Zc.copy(), "PiI": Zc.copy(),
        "Phi": Phi, "Pi_Phi": Zc.copy(),
        "psi": psi, "Pi_psi": Zc.copy(),
        "ptum_params": {"V0": 1.0, "lam": 0.05, "m2": m2,
                        "alpha_p": alpha_p, "Eth": Eth, "lam_psi": lam_psi},
    }


__all__ = ["kerr_complex_ptum_state"]
