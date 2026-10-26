"""Quick AMR long run test - just a few steps to verify pipeline."""
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


def test_amr_pipeline(use_projection=True, n_steps=10):
    """Test the AMR pipeline with a few steps."""
    print(f"Testing AMR pipeline (projection={'ON' if use_projection else 'OFF'})...")

    # Small config for quick test
    N_c = 25
    L_c = 40.0
    N_f = 13
    L_f = 8.0
    dt = 0.125

    h_c = L_c / (N_c - 1)
    h_f = L_f / (N_f - 1)

    # Initial data
    s_c_base = bowen_york_binary(N_c, L_c, m1=1.0, m2=1.0, d=10.0)
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

    # Fine grid mask - use the projection band mask
    band_mask = _band_mask(s_f.phi, h_f, n_layers=4, mask_threshold=0.5)

    H_f_band_history = []
    H_c_max_history = []

    for step in range(n_steps):
        t = step * dt

        # 1. Evolve coarse grid
        s_c = rk4_step_mp_vec_damp(s_c, h_c, dt, eta=1.0, width_cells=3, kappa=0.5)

        # Check for NaN in coarse grid
        if np.any(np.isnan(s_c.phi)):
            print(f"  Coarse grid NaN at step {step}!")
            break

        # 2. Update fine patch position
        c_left, c_right, _ = peak_centroids(s_c, L_c)
        com_x = (c_left[0] + c_right[0]) / 2.0

        # Check for NaN in COM
        if np.isnan(com_x):
            print(f"  COM NaN at step {step}!")
            break

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

        # 4. Apply projection
        if use_projection:
            H_f = hamiltonian_constraint_full(s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
            phi_new, dphi_max = project_phi_band_covariant(
                s_f.phi, s_f.gtilde, h_f, H_f,
                n_layers=4, mask_threshold=0.5
            )
            At_new, l_before, l_after = momentum_band_correct(
                s_f.Atilde, s_f.gtilde, h_f,
                n_layers=4, passes=1
            )

            s_f = BSSNStateMP(
                phi=phi_new,
                gtilde=s_f.gtilde,
                K=s_f.K,
                Atilde=At_new,
                Gamma=s_f.Gamma,
                alpha=s_f.alpha,
                beta=s_f.beta,
                B=s_f.B,
            )
        else:
            dphi_max = 0.0
            l_before = 0.0
            l_after = 0.0

        # Diagnostics
        H_c = hamiltonian_constraint_full(s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c)
        H_f = hamiltonian_constraint_full(s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)

        c_left_f, c_right_f, _ = peak_centroids(s_f, L_f)
        sep_f = pair_separation(c_left_f[0], c_right_f[0])

        H_f_band_max = float(np.max(np.abs(H_f[band_mask]))) if band_mask.any() else 0.0
        H_f_band_history.append(H_f_band_max)
        H_c_max_history.append(float(np.max(np.abs(H_c))))

        if step % 5 == 0:
            print(f"Step {step:2d} t={t:5.2f} | "
                  f"H_c_max={np.max(np.abs(H_c)):.3e} "
                  f"H_f_max={np.max(np.abs(H_f)):.3e} "
                  f"H_f_band={H_f_band_max:.3e} "
                  f"sep={sep_f:.3f} "
                  f"com={com_x:.3f} "
                  f"dphi={dphi_max:.3e} "
                  f"l_before={l_before:.3e} l_after={l_after:.3e}")

        if H_f_band_max > 10.0:
            print(f"  BLOWUP detected at step {step}!")
            break

    return H_f_band_history, H_c_max_history


if __name__ == "__main__":
    # Run with projection
    print("=" * 60)
    hist_proj, hc_proj = test_amr_pipeline(use_projection=True, n_steps=10)

    print("\n" + "=" * 60)
    # Run without projection
    hist_noproj, hc_noproj = test_amr_pipeline(use_projection=False, n_steps=10)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"With projection:    H_f_band_max final = {hist_proj[-1]:.3e}")
    print(f"Without projection: H_f_band_max final = {hist_noproj[-1]:.3e}")
    if hist_proj[-1] > 0:
        print(f"Ratio: {hist_noproj[-1] / hist_proj[-1]:.2f}x")

    print(f"\nCoarse grid max H (proj):    {hc_proj[-1]:.3e}")
    print(f"Coarse grid max H (no proj): {hc_noproj[-1]:.3e}")

    if hist_proj[-1] < hist_noproj[-1]:
        print("[OK] PROJECTION HELPS: constraint violation reduced")
    else:
        print("[FAIL] PROJECTION DID NOT HELP: constraint violation not reduced")