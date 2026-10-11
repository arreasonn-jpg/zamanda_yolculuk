"""Two-grid refined BSSN driver (v6.0-alpha5).

Coupling logic per RK4 step:
1. Evolve coarse (bssn_rk4_matter on coarse)
2. Prolongate coarse -> fine (all fields) using tensor ops
3. Evolve fine (bssn_rk4_matter on fine)
4. Restrict fine -> coarse interior
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_rk4_matter import rk4_step_matter
from .refinement_tensor import prolongate_state, restrict_state_into
from .bssn_grid import d_central


def _state_from_dict(d):
    return BSSNState(
        phi=d["phi"], gtilde=d["gtilde"], K=d["K"],
        Atilde=d["Atilde"], Gamma=d["Gamma"],
        alpha=d["alpha"], beta=d["beta"],
    )


def _copy_state(s):
    return BSSNState(
        phi=s.phi.copy(), gtilde=s.gtilde.copy(), K=s.K.copy(),
        Atilde=s.Atilde.copy(), Gamma=s.Gamma.copy(),
        alpha=s.alpha.copy(), beta=s.beta.copy(),
    )


def _prolongate_field(f_c, N_f, offset):
    """For scalar matter fields."""
    from .refinement import prolongate
    return prolongate(f_c, N_f, offset)


def _restrict_field_into(f_f, out_c, N_c, offset):
    from .refinement import restrict
    ox, oy, oz = offset
    N_f = f_f.shape[0]
    half = N_f // 2
    a, b = ox, ox + half - 1
    sl = (slice(a, b+1), slice(a, b+1), slice(a, b+1))
    out_c[sl] = restrict(f_f, N_c, offset)[sl]


def evolve_refined_bssn(state_c, phi_c, Pi_c,
                        state_f, phi_f, Pi_f,
                        h_c, h_f, offset, dt, n_steps,
                        V_func=None):
    """Two-grid BSSN + scalar evolution.

    state_c/state_f: BSSNState on coarse/fine grids
    phi_c/Pi_c, phi_f/Pi_f: scalar fields on coarse/fine grids
    """
    N_c = state_c.phi.shape[0]
    N_f = state_f.phi.shape[0]

    # work on copies
    s_c = _copy_state(state_c)
    pc, Pic = phi_c.copy(), Pi_c.copy()
    s_f = _copy_state(state_f)
    pf, Pif = phi_f.copy(), Pi_f.copy()

    for step in range(n_steps):
        # 1. Coarse evolution
        s_c, pc, Pic = rk4_step_matter(s_c, pc, Pic, h_c, dt, V_func=V_func)

        # 2. Prolongate coarse -> fine
        fine_dict = prolongate_state(s_c, N_f, offset)
        pf_bc = _prolongate_field(pc, N_f, offset)
        Pif_bc = _prolongate_field(Pic, N_f, offset)

        # 3. Override fine boundary from coarse
        for edge in [0, -1]:
            for key in ("phi", "K", "alpha"):
                fine_dict[key][edge, :, :] = fine_dict[key][edge, :, :]
            pf_bc[edge, :, :] = pf_bc[edge, :, :]
            Pif_bc[edge, :, :] = Pif_bc[edge, :, :]

        # Use prolonged state as fine state
        s_f = _state_from_dict(fine_dict)
        pf = pf_bc.copy()
        Pif = Pif_bc.copy()

        # 4. Fine evolution
        s_f, pf, Pif = rk4_step_matter(s_f, pf, Pif, h_f, dt, V_func=V_func)

        # 5. Restrict fine -> coarse interior
        coarse_dict = {
            "phi": s_c.phi, "gtilde": s_c.gtilde, "K": s_c.K,
            "Atilde": s_c.Atilde, "Gamma": s_c.Gamma,
            "alpha": s_c.alpha, "beta": s_c.beta,
        }
        coarse_dict = restrict_state_into(s_f, coarse_dict, N_c, offset)
        s_c = _state_from_dict(coarse_dict)
        _restrict_field_into(pf, pc, N_c, offset)
        _restrict_field_into(Pif, Pic, N_c, offset)

    return {"state_coarse": s_c, "phi_coarse": pc, "Pi_coarse": Pic,
            "state_fine": s_f, "phi_fine": pf, "Pi_fine": Pif}


__all__ = ["evolve_refined_bssn"]
