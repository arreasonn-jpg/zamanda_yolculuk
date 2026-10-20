"""Refined full stack: Kerr+SC + condensate + PTUM (v6.5-alpha2).

2-grid (coarse+fine) evolution of:
- Kerr+SC metric
- SC condensate (phi_sc, Pi_sc)
- PTUM fields (Phi, Pi_Phi, psi, Pi_psi)

Uses refinement_tensor for all BSSN fields + prolongate/restrict
for the 6 matter fields.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rk4_sc_ptum import rk4_step_sc_ptum
from .amr.refinement_tensor import prolongate_state, restrict_state_into
from .amr.refinement import prolongate, restrict


def _copy_state(s):
    return BSSNState(
        phi=s.phi.copy(), gtilde=s.gtilde.copy(), K=s.K.copy(),
        Atilde=s.Atilde.copy(), Gamma=s.Gamma.copy(),
        alpha=s.alpha.copy(), beta=s.beta.copy(),
    )


def _state_from_dict(d):
    return BSSNState(
        phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
        Atilde=d["Atilde"], Gamma=d["Gamma"],
        alpha=d["alpha"], beta=d["beta"],
    )


def _prolong(f_c, N_f, offset):
    return prolongate(f_c, N_f, offset)


def _restrict_into(f_f, out_c, N_c, offset):
    r = restrict(f_f, N_c, offset)
    valid = ~np.isnan(r)
    out_c[valid] = r[valid]


def evolve_refined_sc_ptum(s_c, psc_c, Psc_c, Ph_c, PPh_c, ps_c, Pps_c,
                            s_f, psc_f, Psc_f, Ph_f, PPh_f, ps_f, Pps_f,
                            h_c, h_f, offset, dt, n_steps,
                            V_sc, ptum_params):
    N_c = s_c.phi.shape[0]
    N_f = s_f.phi.shape[0]

    sc = _copy_state(s_c)
    psc, Psc = psc_c.copy(), Psc_c.copy()
    Ph, PPh = Ph_c.copy(), PPh_c.copy()
    ps, Pps = ps_c.copy(), Pps_c.copy()

    sf = _copy_state(s_f)
    psf, Psf = psc_f.copy(), Psc_f.copy()
    Phf, PPhf = Ph_f.copy(), PPh_f.copy()
    psf_, Ppsf = ps_f.copy(), Pps_f.copy()

    for step in range(n_steps):
        # 1. Coarse RK4
        sc, psc, Psc, Ph, PPh, ps, Pps = rk4_step_sc_ptum(
            sc, psc, Psc, Ph, PPh, ps, Pps, h_c, dt, V_sc, ptum_params)

        # 2. Prolongate coarse -> fine
        fine_dict = prolongate_state(sc, N_f, offset)
        psf = _prolong(psc, N_f, offset)
        Psf = _prolong(Psc, N_f, offset)
        Phf = _prolong(Ph, N_f, offset)
        PPhf = _prolong(PPh, N_f, offset)
        psf_ = _prolong(ps, N_f, offset)
        Ppsf = _prolong(Pps, N_f, offset)

        sf = _state_from_dict(fine_dict)

        # 3. Fine RK4
        sf, psf, Psf, Phf, PPhf, psf_, Ppsf = rk4_step_sc_ptum(
            sf, psf, Psf, Phf, PPhf, psf_, Ppsf,
            h_f, dt, V_sc, ptum_params)

        # 4. Restrict fine -> coarse
        coarse_dict = {
            "phi": sc.phi, "gtilde": sc.gtilde, "K": sc.K,
            "Atilde": sc.Atilde, "Gamma": sc.Gamma,
            "alpha": sc.alpha, "beta": sc.beta,
        }
        coarse_dict = restrict_state_into(sf, coarse_dict, N_c, offset)
        sc = _state_from_dict(coarse_dict)
        _restrict_into(psf, psc, N_c, offset)
        _restrict_into(Psf, Psc, N_c, offset)
        _restrict_into(Phf, Ph, N_c, offset)
        _restrict_into(PPhf, PPh, N_c, offset)
        _restrict_into(psf_, ps, N_c, offset)
        _restrict_into(Ppsf, Pps, N_c, offset)

    return {"state_c": sc, "phi_sc_c": psc, "Pi_sc_c": Psc,
            "Phi_c": Ph, "Pi_Phi_c": PPh, "psi_c": ps, "Pi_psi_c": Pps,
            "state_f": sf, "phi_sc_f": psf, "Pi_sc_f": Psf,
            "Phi_f": Phf, "Pi_Phi_f": PPhf, "psi_f": psf_, "Pi_psi_f": Ppsf}


__all__ = ["evolve_refined_sc_ptum"]
