"""v13: moving fine-patch AMR (window recentering on BH COM).

v12.9-v12.15 verdict chain: fixed-patch AMR loses constraint control
over time; cubic prolongation (v12.11), finer fine (v12.13) and
Gamma-consistent injection (v12.15) are all not the lever. The
remaining physical suspicion (v12.8): the punctures must be properly
resolved, and a FIXED window cannot follow the binary.

Mechanism here: the fine window is recentered (full re-extract +
prolongate from the current coarse state) whenever the BH center of
mass drifts by more than 1 coarse cell from the window center,
snapped to even coarse offsets so the factor-2 prolongation grid
stays aligned.

Injection as the v12.9 baseline (prolongated Gamma; v12.15 showed
Gamma rederivation is not the lever).

Diagnostics (same conventions as the v12.15 baseline run):
- coarse H: fixed centered sub-window mask, H0_c at step 0
- fine H: mask follows the CURRENT fine window (R > 1.5 in window
  coords), normalized with H0_f at step 0 (initial, centered)
- metric of comparison: H_f/H0 @ step 1000 vs baseline 1.0756
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector)
from zt005.amr.amr_bh_track import find_bh_centers
from zt005.amr.amr_moving import target_lo as target_lo_impl
from zt005.amr.amr_moving import fine_grid_mask

N_c, L_c = 17, 20.0
N_f, L_f = 17, 10.0
N_sub = 9
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
DT = 2e-3
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
CENTER_IDX = (N_c - N_sub) // 2


def extract_window(s_c, lo):
    hi = lo + N_sub
    return BSSNStateMP(
        phi=prolongate_scalar(s_c.phi[lo:hi, lo:hi, lo:hi], 2),
        gtilde=prolongate_tensor(s_c.gtilde[lo:hi, lo:hi, lo:hi], 2),
        K=prolongate_scalar(s_c.K[lo:hi, lo:hi, lo:hi], 2),
        Atilde=prolongate_tensor(s_c.Atilde[lo:hi, lo:hi, lo:hi], 2),
        Gamma=prolongate_vector(s_c.Gamma[lo:hi, lo:hi, lo:hi], 2),
        alpha=prolongate_scalar(s_c.alpha[lo:hi, lo:hi, lo:hi], 2),
        beta=prolongate_vector(s_c.beta[lo:hi, lo:hi, lo:hi], 2),
        B=prolongate_vector(s_c.B[lo:hi, lo:hi, lo:hi], 2),
    )


def _band(fine, coarse, n_layers):
    out = fine.copy()
    for k in range(n_layers):
        out[k, ...] = coarse[k, ...]
        out[-1 - k, ...] = coarse[-1 - k, ...]
        out[:, k, ...] = coarse[:, k, ...]
        out[:, -1 - k, ...] = coarse[:, -1 - k, ...]
        out[:, :, k, ...] = coarse[:, :, k, ...]
        out[:, :, -1 - k, ...] = coarse[:, :, -1 - k, ...]
    return out


def inject_full(s_old, s_fresh, n_layers=4):
    return BSSNStateMP(
        phi=_band(s_old.phi, s_fresh.phi, n_layers),
        gtilde=_band(s_old.gtilde, s_fresh.gtilde, n_layers),
        K=_band(s_old.K, s_fresh.K, n_layers),
        Atilde=_band(s_old.Atilde, s_fresh.Atilde, n_layers),
        Gamma=_band(s_old.Gamma, s_fresh.Gamma, n_layers),
        alpha=_band(s_old.alpha, s_fresh.alpha, n_layers),
        beta=_band(s_old.beta, s_fresh.beta, n_layers),
        B=_band(s_old.B, s_fresh.B, n_layers),
    )


def fine_mask():
    return fine_grid_mask(N_f, L_f, 1.5)


def main():
    s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
    lo = target_lo_impl(0.0, h_c, N_c, N_sub)   # initial: centered
    s_f = extract_window(s_c, lo)
    mask_f = fine_mask()

    # fixed coarse diagnostics (centered window, as the v12.15 baseline)
    hi_c = CENTER_IDX + N_sub
    xc = np.linspace(-L_c / 2, L_c / 2, N_c)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    mask_c_sub = (np.sqrt(Xc ** 2 + Yc ** 2 + Zc ** 2)
                  > 1.5)[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c]
    H0_c = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_c.phi[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        s_c.gtilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        s_c.K[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        s_c.Atilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        h_c)[mask_c_sub])))

    # fixed fine normalization: centered window at step 0
    H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

    t_hist, H_c_hist, H_f_hist = [], [], []
    recenter_steps = []
    t0 = time.time()
    for step in range(1, STEPS + 1):
        s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        xl, xr, _ = find_bh_centers(s_c, L_c, N_c, axis=0)
        lo_new = target_lo_impl(0.5 * (xl + xr), h_c, N_c, N_sub)
        if lo_new != lo:
            lo = lo_new
            s_f = extract_window(s_c, lo)
            recenter_steps.append(step)
        s_fresh = extract_window(s_c, lo)
        s_f = inject_full(s_f, s_fresh, n_layers=4)
        s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        if step % 100 == 0 or step == STEPS:
            H_c = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_c.phi[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
                s_c.gtilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                           CENTER_IDX:hi_c],
                s_c.K[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
                s_c.Atilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                           CENTER_IDX:hi_c], h_c)[mask_c_sub])))
            H_f = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
            t_hist.append(step * DT)
            H_c_hist.append(H_c / H0_c)
            H_f_hist.append(H_f / H0_f)
            fin = np.all(np.isfinite(s_c.phi)) \
                and np.all(np.isfinite(s_f.phi))
            print(f"[moving] step={step:4d} lo={lo} "
                  f"H_c/H0={H_c / H0_c:.4f} H_f/H0={H_f / H0_f:.4f} "
                  f"recenters={len(recenter_steps)} "
                  f"fin={fin} t={time.time() - t0:.0f}s", flush=True)
            if not fin:
                print("NON-FINITE"); break

    out = {"steps": STEPS, "recenter_steps": recenter_steps,
           "recenters": len(recenter_steps),
           "H0_fine": H0_f, "H0_coarse_sub": H0_c,
           "t": t_hist, "H_coarse_sub": H_c_hist, "H_fine": H_f_hist}
    base = Path("results/validated")
    np.savez(base / "v13_amr_moving.npz",
             t=np.array(t_hist), H_coarse_sub=np.array(H_c_hist),
             H_fine=np.array(H_f_hist))
    (base / "v13_amr_moving.json").write_text(json.dumps({
        "stage": "v13.0.0", "steps": STEPS,
        "recenter_steps": recenter_steps,
        "recenters": len(recenter_steps),
        "H_fine_final": H_f_hist[-1],
        "baseline_reference_v12_15_step1000": 1.0756,
    }, indent=2))
    print(json.dumps({"H_fine_final": H_f_hist[-1],
                      "H_coarse_final": H_c_hist[-1],
                      "recenters": len(recenter_steps),
                      "recenter_steps": recenter_steps}, indent=2))


if __name__ == "__main__":
    main()
