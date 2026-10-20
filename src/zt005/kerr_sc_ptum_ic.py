"""Kerr + SC metric + SC condensate + PTUM fields (v6.4-alpha1).

Full stack:
- Metric: Kerr-Schild + smoothed M_eff (kerr_sc_ic)
- SC condensate: single real scalar, perturbative amplitude
- PTUM fields: fuel Phi + proxy psi, physical amplitude

Diagnostics compute the total stress-energy (sum of SC + PTUM).
"""
from __future__ import annotations
import numpy as np
from .kerr_sc_ic import kerr_sc_state
from .matter.sc_condensate import sc_condensate_profile, sc_condensate_potential
from .matter.ptum_field import fuel_profile, psi_equilibrium
from .matter.matter_scalar import scalar_sources
from .matter.matter_ptum import ptum_sources


def kerr_sc_ptum_state(N, L, M=1.0, a=0.0, rc=1.0,
                        # SC condensate
                        rho_core=1e19, m_pair_MeV=10.0, g_int=0.1,
                        phi_amp_target=1e-4,
                        # PTUM
                        A_fuel=3.0, sigma_fuel=None,
                        m2=1.0, alpha_p=1.0, Eth=1.0, lam_psi=1.0):
    """Return dict with state + all fields + potentials."""
    s = kerr_sc_state(N, L, M=M, a=a, rc=rc)

    # SC condensate
    phi_sc, Pi_sc, sigma_sc, omega_sc, A_amp = sc_condensate_profile(
        N=N, L=L, M=M, rho_core=rho_core, m_pair_MeV=m_pair_MeV)
    maxv = np.max(np.abs(phi_sc))
    if maxv > 0:
        phi_sc = phi_sc * (phi_amp_target / maxv)
    V_sc = sc_condensate_potential(m_pair_MeV=m_pair_MeV, g_int=g_int)

    # PTUM fields
    if sigma_fuel is None:
        sigma_fuel = L / 5.0
    Phi = fuel_profile(N, L, Phi0=A_fuel, sigma=sigma_fuel)
    psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
    Z = np.zeros((N, N, N))

    return {
        "state": s,
        # SC
        "phi_sc": phi_sc, "Pi_sc": Pi_sc, "V_sc": V_sc,
        # PTUM
        "Phi": Phi, "Pi_Phi": Z.copy(), "psi": psi, "Pi_psi": Z.copy(),
        "V_ptum_params": {"V0": 1.0, "lam": 0.05, "m2": m2,
                          "alpha_p": alpha_p, "Eth": Eth, "lam_psi": lam_psi},
    }


def combined_diagnostics(pack, h):
    """Return total rho, split rho_SC / rho_PTUM, negative cells etc."""
    s = pack["state"]
    # SC
    src_sc = scalar_sources(pack["phi_sc"], pack["Pi_sc"],
                             s.gtilde, s.phi, h, V_func=pack["V_sc"])
    # PTUM
    p = pack["V_ptum_params"]
    src_pt = ptum_sources(pack["Phi"], pack["Pi_Phi"],
                           pack["psi"], pack["Pi_psi"],
                           s.gtilde, s.phi, h,
                           V0=p["V0"], lam=p["lam"], m2=p["m2"],
                           alpha=p["alpha_p"], Eth=p["Eth"],
                           lam_psi=p["lam_psi"])
    rho_tot = src_sc["rho"] + src_pt["rho"]
    return {
        "rho_total_min": float(rho_tot.min()),
        "rho_sc_min": float(src_sc["rho"].min()),
        "rho_ptum_min": float(src_pt["rho"].min()),
        "n_neg_total": int(np.sum(rho_tot < 0)),
        "n_neg_sc": int(np.sum(src_sc["rho"] < 0)),
        "n_neg_ptum": int(np.sum(src_pt["rho"] < 0)),
        "Kmax": float(np.max(np.abs(s.K))),
        "alpha_min": float(s.alpha.min()),
    }


__all__ = ["kerr_sc_ptum_state", "combined_diagnostics"]
