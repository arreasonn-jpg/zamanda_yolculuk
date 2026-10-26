"""Minimal AMR single step test - verify projection works."""
from __future__ import annotations

import numpy as np
from pathlib import Path

from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_rk4_mp_vec_damp import rk4_step_mp_vec_damp
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr.refinement_tensor import prolongate_state
from zt005.amr.amr_interp import inject_boundary
from zt005.amr.amr_moving import target_lo
from zt005.amr.amr_bh_track import peak_centroids, pair_separation
from zt005.amr.amr_projection import (
    project_phi_band_covariant,
    momentum_band_correct,
    _band_mask
)
from zt005.amr.refinement import prolongate


def test_single_step():
    """Test a single AMR step with projection."""
    print("Testing single AMR step with projection...")

    # Very small config for quick test
    N_c = 17
    L_c = 24.0
    N_f = 9
    L_f = 6.0
    dt = 0.125

    h_c = L_c / (N_c - 1)
    h_f = L_f / (N_f - 1)

    # Initial data
    s_c_base = bowen_york_binary(N_c, L_c, m1=1.0, m2=1.0, d=8.0)
    s_c = BSSNStateMP.from_base(s_c_base)

    # Prolongate to fine grid
    offset = ((N_c - N_f) // 2,) * 3
    fine_dict = prolongate_state(s_c, N_f, offset)
    B_f = np.zeros((N_f, N_f, N_f, 3))
    for i in range(3):
        B_f[..., i] = prolongate(s_c.B[..., i], N_f, offset)
    fine_dict["B"] = B_f

    s_f = BSSNStateMP(
        phi=fine_dict["phi"],
        gtilde=fine_dict["gtilde"],
        K=fine_dict["K"],
        Atilde=fine_dict["Atilde"],
        Gamma=fine_dict["Gamma"],
        alpha=fine_dict["alpha"],
        beta=fine_dict["beta"],
        B=fine_dict["B"],
    )

    print(f"Initial H_c max: {np.max(np.abs(hamiltonian_constraint_full(s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c))):.3e}")
    print(f"Initial H_f max: {np.max(np.abs(hamiltonian_constraint_full(s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f))):.3e}")

    # Fine grid mask
    band_mask = _band_mask(s_f.phi, h_f, n_layers=4, mask_threshold=0.5)
    print(f"Band mask cells: {np.sum(band_mask)}")

    # 1. Evolve coarse grid
    s_c = rk4_step_mp_vec_damp(s_c, h_c, dt, eta=1.0, width_cells=3, kappa=0.5)
    print(f"After coarse step: H_c max = {np.max(np.abs(hamiltonian_constraint_full(s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c))):.3e}")

    # 2. Update fine patch position
    c_left, c_right, _ = peak_centroids(s_c, L_c)
    com_x = (c_left[0] + c_right[0]) / 2.0
    lo = target_lo(com_x, h_c, N_c, N_f)
    offset = (lo, lo, lo)

    # Prolongate coarse -> fine at new position
    fine_dict = prolongate_state(s_c, N_f, offset)
    B_f = np.zeros((N_f, N_f, N_f, 3))
    for i in range(3):
        B_f[..., i] = prolongate(s_c.B[..., i], N_f, offset)
    fine_dict["B"] = B_f

    # Inject boundary
    for key in ["phi", "K", "alpha", "gtilde", "Atilde", "Gamma", "beta"]:
        fine_dict[key] = inject_boundary(fine_dict[key], fine_dict[key], n_layers=4)

    s_f = BSSNStateMP(
        phi=fine_dict["phi"],
        gtilde=fine_dict["gtilde"],
        K=fine_dict["K"],
        Atilde=fine_dict["Atilde"],
        Gamma=fine_dict["Gamma"],
        alpha=fine_dict["alpha"],
        beta=fine_dict["beta"],
        B=fine_dict["B"],
    )

    # 3. Evolve fine grid
    s_f = rk4_step_mp_vec_damp(s_f, h_f, dt, eta=1.0, width_cells=3, kappa=0.5)
    print(f"After fine step: H_f max = {np.max(np.abs(hamiltonian_constraint_full(s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f))):.3e}")

    # 4. Apply projection
    H_f = hamiltonian_constraint_full(s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
    H_f_band_before = float(np.max(np.abs(H_f[band_mask]))) if band_mask.any() else 0.0
    print(f"H_f band before projection: {H_f_band_before:.3e}")

    phi_new, dphi_max = project_phi_band_covariant(
        s_f.phi, s_f.gtilde, h_f, H_f,
        n_layers=4, mask_threshold=0.5
    )
    At_new, l_before, l_after = momentum_band_correct(
        s_f.Atilde, s_f.gtilde, h_f,
        n_layers=4, passes=1
    )

    s_f_proj = BSSNStateMP(
        phi=phi_new,
        gtilde=s_f.gtilde,
        K=s_f.K,
        Atilde=At_new,
        Gamma=s_f.Gamma,
        alpha=s_f.alpha,
        beta=s_f.beta,
        B=s_f.B,
    )

    H_f_proj = hamiltonian_constraint_full(s_f_proj.phi, s_f_proj.gtilde, s_f_proj.K, s_f_proj.Atilde, h_f)
    H_f_band_after = float(np.max(np.abs(H_f_proj[band_mask]))) if band_mask.any() else 0.0
    print(f"H_f band after projection: {H_f_band_after:.3e}")
    print(f"dphi_max: {dphi_max:.3e}")
    print(f"l_before: {l_before:.3e}, l_after: {l_after:.3e}")

    if H_f_band_after < H_f_band_before:
        print("[OK] PROJECTION REDUCES BAND CONSTRAINT VIOLATION")
    else:
        print("[FAIL] PROJECTION DID NOT REDUCE BAND CONSTRAINT VIOLATION")

    return H_f_band_before, H_f_band_after


if __name__ == "__main__":
    test_single_step()