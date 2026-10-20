"""Two-grid refined BSSN + PTUM driver (v6.0-alpha6).

Coarse: N_c, h_c
Fine:   N_f, h_f = h_c/2 (centered)
Matter: 4 real fields (Phi, Pi_Phi, psi, Pi_psi)
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rk4_ptum import rk4_step_ptum
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


def _prolong_scalar(f_c, N_f, offset):
    return prolongate(f_c, N_f, offset)


def _restrict_into_scalar(f_f, out_c, N_c, offset):
    r = restrict(f_f, N_c, offset)
    valid = ~np.isnan(r)
    out_c[valid] = r[valid]


def evolve_refined_ptum(s_c, Phi_c, PiPhi_c, psi_c, Pipsi_c,
                        s_f, Phi_f, PiPhi_f, psi_f, Pipsi_f,
                        h_c, h_f, offset, dt, n_steps,
                        V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0,
                        lam_psi=1.0):
    N_c = s_c.phi.shape[0]
    N_f = s_f.phi.shape[0]

    sc = _copy_state(s_c)
    Phc, PPhc = Phi_c.copy(), PiPhi_c.copy()
    psc, Ppsc = psi_c.copy(), Pipsi_c.copy()

    sf = _copy_state(s_f)
    Phf, PPhf = Phi_f.copy(), PiPhi_f.copy()
    psf, Ppsf = psi_f.copy(), Pipsi_f.copy()

    for step in range(n_steps):
        # 1. Coarse PTUM RK4
        sc, Phc, PPhc, psc, Ppsc = rk4_step_ptum(
            sc, Phc, PPhc, psc, Ppsc, h_c, dt,
            V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth,
            lam_psi=lam_psi)

        # 2. Prolongate coarse -> fine
        fine_dict = prolongate_state(sc, N_f, offset)
        Phf = _prolong_scalar(Phc, N_f, offset)
        PPhf = _prolong_scalar(PPhc, N_f, offset)
        psf = _prolong_scalar(psc, N_f, offset)
        Ppsf = _prolong_scalar(Ppsc, N_f, offset)

        sf = _state_from_dict(fine_dict)

        # 3. Fine PTUM RK4
        sf, Phf, PPhf, psf, Ppsf = rk4_step_ptum(
            sf, Phf, PPhf, psf, Ppsf, h_f, dt,
            V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth,
            lam_psi=lam_psi)

        # 4. Restrict fine -> coarse interior
        coarse_dict = {
            "phi": sc.phi, "gtilde": sc.gtilde, "K": sc.K,
            "Atilde": sc.Atilde, "Gamma": sc.Gamma,
            "alpha": sc.alpha, "beta": sc.beta,
        }
        coarse_dict = restrict_state_into(sf, coarse_dict, N_c, offset)
        sc = _state_from_dict(coarse_dict)
        _restrict_into_scalar(Phf, Phc, N_c, offset)
        _restrict_into_scalar(PPhf, PPhc, N_c, offset)
        _restrict_into_scalar(psf, psc, N_c, offset)
        _restrict_into_scalar(Ppsf, Ppsc, N_c, offset)

    return {"state_c": sc, "Phi_c": Phc, "PiPhi_c": PPhc,
            "psi_c": psc, "Pipsi_c": Ppsc,
            "state_f": sf, "Phi_f": Phf, "PiPhi_f": PPhf,
            "psi_f": psf, "Pipsi_f": Ppsf}


__all__ = ["evolve_refined_ptum"]
