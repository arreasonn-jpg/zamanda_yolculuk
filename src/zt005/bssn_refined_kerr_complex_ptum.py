"""Refined Kerr + complex scalar + PTUM driver (v7.1-alpha2)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rk4_complex_ptum import rk4_step_complex_ptum
from .amr.refinement_tensor import prolongate_state, restrict_state_into
from .amr.refinement import prolongate
from .matter.ptum_field import psi_equilibrium_smooth


def _copy_state(s):
    return BSSNState(phi=s.phi.copy(), gtilde=s.gtilde.copy(), K=s.K.copy(),
                     Atilde=s.Atilde.copy(), Gamma=s.Gamma.copy(),
                     alpha=s.alpha.copy(), beta=s.beta.copy())


def _state_from_dict(d):
    return BSSNState(phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
                     Atilde=d["Atilde"], Gamma=d["Gamma"],
                     alpha=d["alpha"], beta=d["beta"])


def evolve_refined_kerr_complex_ptum(
        s_c, PRc, PIc, PiRc, PiIc, Phc, PPhc, psc, Ppsc,
        s_f, PRf, PIf, PiRf, PiIf, Phf, PPhf, psf, Ppsf,
        h_c, h_f, offset, dt, n_steps, ptum_params,
        psi_from_equilibrium=True, dV_func=None):
    N_c = s_c.phi.shape[0]
    N_f = s_f.phi.shape[0]

    sc = _copy_state(s_c)
    PR, PI, PiR, PiI = PRc.copy(), PIc.copy(), PiRc.copy(), PiIc.copy()
    Ph, PPh, ps, Pps = Phc.copy(), PPhc.copy(), psc.copy(), Ppsc.copy()

    sf = _copy_state(s_f)
    PRf_, PIf_, PiRf_, PiIf_ = PRf.copy(), PIf.copy(), PiRf.copy(), PiIf.copy()
    Phf_, PPhf_, psf_, Ppsf_ = Phf.copy(), PPhf.copy(), psf.copy(), Ppsf.copy()

    m2 = ptum_params.get("m2", 1.0)
    alpha_p = ptum_params.get("alpha_p", 1.0)
    Eth = ptum_params.get("Eth", 1.0)
    lam_psi = ptum_params.get("lam_psi", 1.0)

    for step in range(n_steps):
        sc, PR, PI, PiR, PiI, Ph, PPh, ps, Pps = rk4_step_complex_ptum(
            sc, PR, PI, PiR, PiI, Ph, PPh, ps, Pps,
            h_c, dt, ptum_params, dV_func)

        fine_dict = prolongate_state(sc, N_f, offset)
        PRf_ = prolongate(PR, N_f, offset)
        PIf_ = prolongate(PI, N_f, offset)
        PiRf_ = prolongate(PiR, N_f, offset)
        PiIf_ = prolongate(PiI, N_f, offset)
        Phf_ = prolongate(Ph, N_f, offset)
        PPhf_ = prolongate(PPh, N_f, offset)
        if psi_from_equilibrium:
            psf_ = psi_equilibrium_smooth(Phf_, m2=m2, alpha=alpha_p,
                                            Eth=Eth, lam_psi=lam_psi)
            Ppsf_ = np.zeros_like(Ppsf_)
        else:
            psf_ = prolongate(ps, N_f, offset)
            Ppsf_ = prolongate(Pps, N_f, offset)
        sf = _state_from_dict(fine_dict)

        sf, PRf_, PIf_, PiRf_, PiIf_, Phf_, PPhf_, psf_, Ppsf_ = \
            rk4_step_complex_ptum(sf, PRf_, PIf_, PiRf_, PiIf_,
                                   Phf_, PPhf_, psf_, Ppsf_,
                                   h_f, dt, ptum_params, dV_func)

        coarse_dict = {"phi": sc.phi, "gtilde": sc.gtilde, "K": sc.K,
                       "Atilde": sc.Atilde, "Gamma": sc.Gamma,
                       "alpha": sc.alpha, "beta": sc.beta}
        coarse_dict = restrict_state_into(sf, coarse_dict, N_c, offset)
        sc = _state_from_dict(coarse_dict)

    return {"state_c": sc, "PR_c": PR, "PI_c": PI, "PiR_c": PiR, "PiI_c": PiI,
            "Ph_c": Ph, "PPh_c": PPh, "ps_c": ps, "Pps_c": Pps,
            "state_f": sf, "PR_f": PRf_, "PI_f": PIf_,
            "PiR_f": PiRf_, "PiI_f": PiIf_,
            "Ph_f": Phf_, "PPh_f": PPhf_, "ps_f": psf_, "Pps_f": Ppsf_}


__all__ = ["evolve_refined_kerr_complex_ptum"]
