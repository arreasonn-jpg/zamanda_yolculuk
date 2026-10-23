"""v13.4: combined band projection (covariant phi + momentum).

Same fixed-patch 2x harness as v12.9/v12.15/v13.2 (N_c = 17, window 9,
factor 2 -> N_f = 17, L = 20/10, dt = 2e-3, 4 injected layers, 1000
steps). Per step:
1. inject, then the covariant phi-band solve as in v13.3;
2. momentum_band_correct (v13.4): vector-Poisson correction of the
   Atilde divergence inside the (face-skipped) band, three passes of
   the documented fixed-point iteration.
Metrics: H_f/H0_raw + H_f/H0_iso as in v13.3, plus the momentum
metric M/i0 = max |lambda| divergence over the inner band, normalized
to the uncorrected step-0 value, and the per-step correction ratio.

References:
- baseline injection (no projection): H_f/H0_raw @ 1000 = 1.0756
- v13.2 flat projection:              H_f/H0_raw @ 1000 = 0.9319
- v13.3 covariant projection:         H_f/H0_raw @ 1000 = 0.9382
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
from zt005.amr.amr_projection import (
    _geometric_band, atilde_up_divergence, momentum_band_correct,
    project_phi_band_covariant)

N_c, L_c = 17, 20.0
N_sub = 9
factor = 2
N_f = (N_sub - 1) * factor + 1
L_f = 10.0
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
DT = 2e-3
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
LO, HI = (N_c - N_sub) // 2, (N_c + N_sub) // 2


def extract_window(s_c):
    return BSSNStateMP(
        phi=prolongate_scalar(s_c.phi[LO:HI, LO:HI, LO:HI], factor),
        gtilde=prolongate_tensor(s_c.gtilde[LO:HI, LO:HI, LO:HI], factor),
        K=prolongate_scalar(s_c.K[LO:HI, LO:HI, LO:HI], factor),
        Atilde=prolongate_tensor(s_c.Atilde[LO:HI, LO:HI, LO:HI], factor),
        Gamma=prolongate_vector(s_c.Gamma[LO:HI, LO:HI, LO:HI], factor),
        alpha=prolongate_scalar(s_c.alpha[LO:HI, LO:HI, LO:HI], factor),
        beta=prolongate_vector(s_c.beta[LO:HI, LO:HI, LO:HI], factor),
        B=prolongate_vector(s_c.B[LO:HI, LO:HI, LO:HI], factor),
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


def inject_full(s_old, s_new_data, n_layers=4):
    return BSSNStateMP(
        phi=_band(s_old.phi, s_new_data.phi, n_layers),
        gtilde=_band(s_old.gtilde, s_new_data.gtilde, n_layers),
        K=_band(s_old.K, s_new_data.K, n_layers),
        Atilde=_band(s_old.Atilde, s_new_data.Atilde, n_layers),
        Gamma=_band(s_old.Gamma, s_new_data.Gamma, n_layers),
        alpha=_band(s_old.alpha, s_new_data.alpha, n_layers),
        beta=_band(s_old.beta, s_new_data.beta, n_layers),
        B=_band(s_old.B, s_new_data.B, n_layers),
    )


def main():
    s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
    s_f = extract_window(s_c)
    xc = np.linspace(-L_c / 2, L_c / 2, N_c)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    R_c = np.sqrt(Xc ** 2 + Yc ** 2 + Zc ** 2)
    mask_c = np.zeros((N_c, N_c, N_c), dtype=bool)
    mask_c[LO:HI, LO:HI, LO:HI] = \
        R_c[LO:HI, LO:HI, LO:HI] > 1.5
    H0_c_full = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c)[mask_c])))
    H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f))))
    # iso-reference: H0 after one covariant projection on the initial
    # data (measurement only; the state keeps the raw prolongation)
    H_res0 = hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
    phi_proj0, _ = project_phi_band_covariant(
        s_f.phi, s_f.gtilde, h_f, H_res0, n_layers=4,
        mask_threshold=0.5)
    H0_f_proj = float(np.mean(np.abs(hamiltonian_constraint_full(
        phi_proj0, s_f.gtilde, s_f.K, s_f.Atilde, h_f))))

    dmax_hist, dvmax_hist, ratio_hist = [], [], []
    t_hist, H_c_hist, H_f_hist, H_fp_hist = [], [], [], []
    M_hist, Mb_hist = [], []
    xf = np.linspace(-L_f / 2, L_f / 2, N_f)
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    mask_f = np.sqrt(Xf ** 2 + Yf ** 2 + Zf ** 2) > 1.5
    # momentum metric references: raw divergence at step 0, global
    # and inner-band-restricted (the band is what the correction can
    # touch; the global max may sit outside it, e.g. near punctures)
    band_m = _geometric_band(N_f, 4)
    M0 = float(np.max(np.abs(atilde_up_divergence(
        s_f.Atilde, s_f.gtilde, h_f))))
    lb0 = atilde_up_divergence(s_f.Atilde, s_f.gtilde, h_f)
    M0_band = float(np.max(np.abs(lb0[band_m])))
    t0 = time.time()
    for step in range(1, STEPS + 1):
        s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        fresh = extract_window(s_c)
        s_f = inject_full(s_f, fresh, n_layers=4)
        # --- v13.3: covariant projection of the phi band ---
        H_res = hamiltonian_constraint_full(
            s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
        phi_new, dmax = project_phi_band_covariant(
            s_f.phi, s_f.gtilde, h_f, H_res, n_layers=4,
            mask_threshold=0.5)
        # --- v13.4: momentum (Atilde divergence) band correction ---
        At_new, l_before, l_after = momentum_band_correct(
            s_f.Atilde, s_f.gtilde, h_f, n_layers=4, sweeps=20,
            omega=0.8, passes=3)
        s_f_corr = BSSNStateMP(
            phi=phi_new, gtilde=s_f.gtilde, K=s_f.K,
            Atilde=At_new, Gamma=s_f.Gamma, alpha=s_f.alpha,
            beta=s_f.beta, B=s_f.B)
        s_f = rk4_step_mp_vec_z4c(s_f_corr, h=h_f, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        dmax_hist.append(dmax)
        dvmax_hist.append(l_before)
        ratio_hist.append(l_after / l_before if l_before > 0 else 1.0)
        if step % 50 == 0 or step == STEPS:
            H_c = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde,
                h_c)[mask_c])))
            H_f = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde,
                h_f)[mask_f])))
            M = float(np.max(np.abs(atilde_up_divergence(
                s_f.Atilde, s_f.gtilde, h_f))))
            l_d = atilde_up_divergence(s_f.Atilde, s_f.gtilde, h_f)
            Mb = float(np.max(np.abs(l_d[band_m])))
            t_hist.append(step * DT)
            H_c_hist.append(H_c / H0_c_full)
            H_f_hist.append(H_f / H0_f)
            H_fp_hist.append(H_f / H0_f_proj)
            M_hist.append(M / M0)
            Mb_hist.append(Mb / M0_band)
            fin = np.all(np.isfinite(s_c.phi)) and \
                np.all(np.isfinite(s_f.phi)) and \
                np.all(np.isfinite(s_f.Atilde))
            print(f"[v13.4 mom] step={step:4d} H_c={H_c/H0_c_full:.4f} "
                  f"H_f={H_f/H0_f:.4f} H_fp={H_f/H0_f_proj:.4f} "
                  f"M/M0={M/M0:.4f} Mb/Mb0={Mb/M0_band:.4f} "
                  f"ratio={ratio_hist[-1]:.3f} "
                  f"dphi={dmax:.1e} fin={fin} t={time.time()-t0:.0f}s",
                  flush=True)
            if not fin:
                print("NON-FINITE")
                break
            if step % 250 == 0:
                base_ckpt = Path("results/validated")
                np.savez(base_ckpt / "v13_4_projection_momentum.npz",
                         t=np.array(t_hist),
                         H_coarse=np.array(H_c_hist),
                         H_fine=np.array(H_f_hist),
                         H_fine_isoref=np.array(H_fp_hist),
                         M_over_M0=np.array(M_hist),
                         Mb_over_Mb0=np.array(Mb_hist),
                         dphi_max=np.array(dmax_hist),
                         dv_max=np.array(dvmax_hist),
                         ratio=np.array(ratio_hist))
                print(f"[checkpoint] v13.4 step={step}", flush=True)

    base = Path("results/validated")
    np.savez(base / "v13_4_projection_momentum.npz",
             t=np.array(t_hist), H_coarse=np.array(H_c_hist),
             H_fine=np.array(H_f_hist),
             H_fine_isoref=np.array(H_fp_hist),
             M_over_M0=np.array(M_hist),
             Mb_over_Mb0=np.array(Mb_hist),
             dphi_max=np.array(dmax_hist),
             dv_max=np.array(dvmax_hist),
             ratio=np.array(ratio_hist))
    (base / "v13_4_projection_momentum.json").write_text(json.dumps({
        "stage": "v13.4",
        "steps": STEPS,
        "H0_fine_raw": H0_f,
        "H0_fine_isoref": H0_f_proj,
        "M0_fine": M0,
        "M0_fine_band": M0_band,
        "H_fine_final_per_raw": H_f_hist[-1],
        "H_fine_final_per_isoref": H_fp_hist[-1],
        "M_final_per_M0": M_hist[-1] if M_hist else None,
        "Mband_final_per_Mb0": Mb_hist[-1] if Mb_hist else None,
        "baseline_reference_step1000": 1.0756,
        "v13_3_reference_step1000": 0.9382,
        "dphi_max_last": dmax_hist[-1] if dmax_hist else None,
        "dvl_ratio_mean": float(np.mean(ratio_hist)),
        "reproduce": "python scripts/v13_4_projection_momentum.py 1000",
    }, indent=2))


if __name__ == "__main__":
    main()
