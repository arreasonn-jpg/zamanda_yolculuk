"""Refined complex scalar + PTUM driver (v6.6-alpha3).

2-grid (coarse + fine) evolution of:
- BSSN metric
- Complex scalar (Psi_R, Psi_I, Pi_R, Pi_I)
- PTUM fields (Phi, Pi_Phi, psi, Pi_psi)
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rk4_complex_ptum import rk4_step_complex_ptum
from .refinement_tensor import prolongate_state, restrict_state_into
from .refinement import prolongate, restrict


def _copy_state(s):
    return BSSNState(phi=s.phi.copy(), gtilde=s.gtilde.copy(), K=s.K.copy(),
                     Atilde=s.Atilde.copy(), Gamma=s.Gamma.copy(),
                     alpha=s.alpha.copy(), beta=s.beta.copy())


def _state_from_dict(d):
    return BSSNState(phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
                     Atilde=d["Atilde"], Gamma=d["Gamma"],
                     alpha=d["alpha"], beta=d["beta"])


def _prolong(f, N_f, offset):
    return prolongate(f, N_f, offset)


def _restrict_into(f_f, out_c, N_c, offset):
    r = restrict(f_f, N_c, offset)
    valid = ~np.isnan(r)
    out_c[valid] = r[valid]


def evolve_refined_complex_ptum(
        s_c, PRc, PIc, PiRc, PiIc, Phc, PPhc, psc, Ppsc,
        s_f, PRf, PIf, PiRf, PiIf, Phf, PPhf, psf, Ppsf,
        h_c, h_f, offset, dt, n_steps, ptum_params, dV_func=None):
    N_c = s_c.phi.shape[0]
    N_f = s_f.phi.shape[0]

    sc = _copy_state(s_c)
    PR, PI, PiR, PiI = PRc.copy(), PIc.copy(), PiRc.copy(), PiIc.copy()
    Ph, PPh, ps, Pps = Phc.copy(), PPhc.copy(), psc.copy(), Ppsc.copy()

    sf = _copy_state(s_f)
    PRf_, PIf_, PiRf_, PiIf_ = PRf.copy(), PIf.copy(), PiRf.copy(), PiIf.copy()
    Phf_, PPhf_, psf_, Ppsf_ = Phf.copy(), PPhf.copy(), psf.copy(), Ppsf.copy()

    for step in range(n_steps):
        # 1. Coarse RK4
        sc, PR, PI, PiR, PiI, Ph, PPh, ps, Pps = rk4_step_complex_ptum(
            sc, PR, PI, PiR, PiI, Ph, PPh, ps, Pps,
            h_c, dt, ptum_params, dV_func)

        # 2. Prolongate coarse -> fine
        fine_dict = prolongate_state(sc, N_f, offset)
        PRf_ = _prolong(PR, N_f, offset)
        PIf_ = _prolong(PI, N_f, offset)
        PiRf_ = _prolong(PiR, N_f, offset)
        PiIf_ = _prolong(PiI, N_f, offset)
        Phf_ = _prolong(Ph, N_f, offset)
        PPhf_ = _prolong(PPh, N_f, offset)
        psf_ = _prolong(ps, N_f, offset)
        Ppsf_ = _prolong(Pps, N_f, offset)
        sf = _state_from_dict(fine_dict)

        # 3. Fine RK4
        sf, PRf_, PIf_, PiRf_, PiIf_, Phf_, PPhf_, psf_, Ppsf_ = \
            rk4_step_complex_ptum(sf, PRf_, PIf_, PiRf_, PiIf_,
                                   Phf_, PPhf_, psf_, Ppsf_,
                                   h_f, dt, ptum_params, dV_func)

        # 4. Restrict fine -> coarse
        coarse_dict = {"phi": sc.phi, "gtilde": sc.gtilde, "K": sc.K,
                       "Atilde": sc.Atilde, "Gamma": sc.Gamma,
                       "alpha": sc.alpha, "beta": sc.beta}
        coarse_dict = restrict_state_into(sf, coarse_dict, N_c, offset)
        sc = _state_from_dict(coarse_dict)
        _restrict_into(PRf_, PR, N_c, offset)
        _restrict_into(PIf_, PI, N_c, offset)
        _restrict_into(PiRf_, PiR, N_c, offset)
        _restrict_into(PiIf_, PiI, N_c, offset)
        _restrict_into(Phf_, Ph, N_c, offset)
        _restrict_into(PPhf_, PPh, N_c, offset)
        _restrict_into(psf_, ps, N_c, offset)
        _restrict_into(Ppsf_, Pps, N_c, offset)

    return {"state_c": sc, "PR_c": PR, "PI_c": PI, "PiR_c": PiR, "PiI_c": PiI,
            "Ph_c": Ph, "PPh_c": PPh, "ps_c": ps, "Pps_c": Pps,
            "state_f": sf, "PR_f": PRf_, "PI_f": PIf_,
            "PiR_f": PiRf_, "PiI_f": PiIf_,
            "Ph_f": Phf_, "PPh_f": PPhf_, "ps_f": psf_, "Pps_f": Ppsf_}


__all__ = ["evolve_refined_complex_ptum"]
