"""v14.1: projection-on LONG evolution of the binary (the marathon).

v14.0 stage-0 verdict: the unprojected moving-patch baseline blows
up at t = 5.6 IN THE AMR INJECTION BAND (fine band max |H|: 0.44 at
t=1 -> 2.32 at t=5, e-fold ~ 0.4 M), while the coarse grid stays
healthy and the orbit-gauge signal is real (omega +0.13..0.15
rad/M near the Newtonian reference 0.125). The full-inspiral
critical path defined there is THIS run: the exact v13.4 projection
protocol (covariant phi-band solve + 3-pass momentum Vector-Poisson
correction) applied every step on the FIXED 2x-patch harness
(no recentering), run LONG (target t = 10, i.e. x4 the validated
window and 0.8x past the unprojected blow-up time).

Instrumentation (from v14.0):
- two puncture locators on the fine window (peak_centroids ring
  meter + alpha_pair basin meter; window offset -L_f/2 + coarse_x[LO]);
- per-500-step H maps (coarse sub-window + fine window);
- full state dump on NON-FINITE (post-mortem support);
- v13.4 metrics: H_f/H0 raw + iso, M/M0, Mb/Mb0, per-step ratio.

Verdict metric: does the projected band stay SUB-CRITICAL to t = 10?
If yes, a 10000-step run (t = 20) becomes the marathon proper.
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
from zt005.amr.amr_bh_track import peak_centroids, alpha_pair

N_c, L_c = 17, 20.0
N_sub = 9
factor = 2
N_f = (N_sub - 1) * factor + 1
L_f = 10.0
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
DT = 2e-3
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
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


def dump_state(path, s_c, s_f, step):
    np.savez(path,
             step=step,
             c_phi=s_c.phi, c_gtilde=s_c.gtilde, c_K=s_c.K,
             c_Atilde=s_c.Atilde, c_Gamma=s_c.Gamma,
             c_alpha=s_c.alpha, c_beta=s_c.beta,
             f_phi=s_f.phi, f_gtilde=s_f.gtilde, f_K=s_f.K,
             f_Atilde=s_f.Atilde, f_Gamma=s_f.Gamma,
             f_alpha=s_f.alpha, f_beta=s_f.beta)


def main():
    s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5,
                                                  4.0))
    s_f = extract_window(s_c)
    xc = np.linspace(-L_c / 2, L_c / 2, N_c)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    R_c = np.sqrt(Xc ** 2 + Yc ** 2 + Zc ** 2)
    mask_c = np.zeros((N_c, N_c, N_c), dtype=bool)
    mask_c[LO:HI, LO:HI, LO:HI] = R_c[LO:HI, LO:HI, LO:HI] > 1.5
    H0_c_full = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c)[mask_c])))
    H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f))))
    H_res0 = hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
    phi_proj0, _ = project_phi_band_covariant(
        s_f.phi, s_f.gtilde, h_f, H_res0, n_layers=4,
        mask_threshold=0.5)
    H0_f_proj = float(np.mean(np.abs(hamiltonian_constraint_full(
        phi_proj0, s_f.gtilde, s_f.K, s_f.Atilde, h_f))))

    xf = np.linspace(-L_f / 2, L_f / 2, N_f)
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    mask_f = np.sqrt(Xf ** 2 + Yf ** 2 + Zf ** 2) > 1.5
    band_m = _geometric_band(N_f, 4)
    M0 = float(np.max(np.abs(atilde_up_divergence(
        s_f.Atilde, s_f.gtilde, h_f))))
    lb0 = atilde_up_divergence(s_f.Atilde, s_f.gtilde, h_f)
    M0_band = float(np.max(np.abs(lb0[band_m])))

    dmax_hist, dvmax_hist, ratio_hist = [], [], []
    t_hist, H_c_hist, H_f_hist, H_fp_hist = [], [], [], []
    M_hist, Mb_hist = [], []
    # v14 instrument series
    cl_hist, cr_hist = [], []
    d_hist, phase_hist, omega_hist = [], [], []
    d_al_hist, amin_hist = [], []
    phase_prev = None
    phase_tot = 0.0
    off = xc[LO] + L_f / 2.0

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
        s_f = rk4_step_mp_vec_z4c(s_f_corr, h=h_f, dt=DT,
                                  width_cells=3, kappa=0.0, eta=1.0)
        dmax_hist.append(dmax)
        dvmax_hist.append(l_before)
        ratio_hist.append(l_after / l_before if l_before > 0 else 1.0)
        if step % 500 == 0:
            H_c_map = hamiltonian_constraint_full(
                s_c.phi, s_c.gtilde, s_c.K, s_c.Atilde, h_c)
            H_f_map = hamiltonian_constraint_full(
                s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)
            np.savez(Path(f"results/validated/"
                          f"v14_1_hmaps_step{step:05d}.npz"),
                     H_coarse_map=H_c_map, H_fine_map=H_f_map,
                     mask_c=mask_c, mask_f=mask_f)
        if step % 50 == 0 or step == STEPS:
            t = step * DT
            (cl_a, cr_a, _) = peak_centroids(s_f, L_f)
            (cl_g, cr_g, depths_g) = alpha_pair(s_f, L_f)
            xl, yl = cl_a[0] + off, cl_a[1] + off
            xr, yr = cr_a[0] + off, cr_a[1] + off
            d = float(np.hypot(xr - xl, yr - yl))
            ang = float(np.arctan2(yr - yl, xr - xl))
            if phase_prev is None:
                dp = 0.0
            else:
                dp = ang - phase_prev
                while dp > np.pi:
                    dp -= 2.0 * np.pi
                while dp < -np.pi:
                    dp += 2.0 * np.pi
                phase_tot += dp
                omega_hist.append(dp / (t - t_hist[-1]))
            phase_prev = ang
            t_hist.append(t)
            phase_hist.append(phase_tot)
            d_hist.append(d)
            d_al_hist.append(float(np.hypot(
                cr_g[0] - cl_g[0], cr_g[1] - cl_g[1])))
            amin_hist.append((float(depths_g[0]),
                              float(depths_g[1])))
            cl_hist.append((xl, yl))
            cr_hist.append((xr, yr))
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
            H_c_hist.append(H_c / H0_c_full)
            H_f_hist.append(H_f / H0_f)
            H_fp_hist.append(H_f / H0_f_proj)
            M_hist.append(M / M0)
            Mb_hist.append(Mb / M0_band)
            om_last = omega_hist[-1] if omega_hist else float("nan")
            fin = np.all(np.isfinite(s_c.phi)) and \
                np.all(np.isfinite(s_f.phi)) and \
                np.all(np.isfinite(s_f.Atilde)) and \
                np.isfinite(d)
            print(f"[v14.1 proj] step={step:4d} t={t:5.2f} "
                  f"H_c={H_c / H0_c_full:.4f} "
                  f"H_f={H_f / H0_f:.4f} "
                  f"H_fp={H_f / H0_f_proj:.4f} "
                  f"M/M0={M / M0:.4f} Mb/Mb0={Mb / M0_band:.4f} "
                  f"ratio={ratio_hist[-1]:.3f} d={d:.3f} "
                  f"phase={phase_tot:+.4f} omega={om_last:+.4f} "
                  f"d_al={d_al_hist[-1]:.3f} "
                  f"a_min={min(depths_g):.3f} fin={fin} "
                  f"WALL={time.time() - t0:.0f}s", flush=True)
            if not fin:
                print("NON-FINITE: dumping state", flush=True)
                dump_state(Path("results/validated/"
                                f"v14_1_state_step{step:05d}.npz"),
                           s_c, s_f, step)
                break
            if step % 500 == 0:
                base_ckpt = Path("results/validated")
                np.savez(base_ckpt / "v14_1_long_projection.npz",
                         t=np.array(t_hist),
                         H_coarse=np.array(H_c_hist),
                         H_fine=np.array(H_f_hist),
                         H_fine_isoref=np.array(H_fp_hist),
                         M_over_M0=np.array(M_hist),
                         Mb_over_Mb0=np.array(Mb_hist),
                         dphi_max=np.array(dmax_hist),
                         dv_max=np.array(dvmax_hist),
                         ratio=np.array(ratio_hist),
                         d=np.array(d_hist),
                         phase=np.array(phase_hist),
                         omega=np.array(omega_hist),
                         d_alpha=np.array(d_al_hist),
                         alpha_min=np.array(amin_hist))
                print(f"[checkpoint] v14.1 step={step}", flush=True)

    base = Path("results/validated")
    cl_arr = np.array(cl_hist)
    cr_arr = np.array(cr_hist)
    np.savez(base / "v14_1_long_projection.npz",
             t=np.array(t_hist),
             H_coarse=np.array(H_c_hist),
             H_fine=np.array(H_f_hist),
             H_fine_isoref=np.array(H_fp_hist),
             M_over_M0=np.array(M_hist),
             Mb_over_Mb0=np.array(Mb_hist),
             dphi_max=np.array(dmax_hist),
             dv_max=np.array(dvmax_hist),
             ratio=np.array(ratio_hist),
             d=np.array(d_hist),
             phase=np.array(phase_hist),
             omega=np.array(omega_hist),
             d_alpha=np.array(d_al_hist),
             alpha_min=np.array(amin_hist),
             cl_x=cl_arr[:, 0], cl_y=cl_arr[:, 1],
             cr_x=cr_arr[:, 0], cr_y=cr_arr[:, 1])
    mean_om = float(np.mean(omega_hist)) if omega_hist else 0.0
    (base / "v14_1_long_projection.json").write_text(json.dumps({
        "stage": "v14.1",
        "steps": STEPS,
        "t_final": float(t_hist[-1]),
        "H_fine_final_per_raw": float(H_f_hist[-1]),
        "H_fine_final_per_isoref": float(H_fp_hist[-1]),
        "M_final_per_M0": float(M_hist[-1]),
        "Mband_final_per_Mb0": float(Mb_hist[-1]),
        "d_ratio": float(d_hist[-1] / 4.0),
        "orbits": phase_tot / (2.0 * np.pi),
        "omega_measured_mean": mean_om,
        "omega_newtonian_ref": 0.125,
        "ratio_mean": float(np.mean(ratio_hist)),
        "v14_0_unprojected_blowup_t": 5.6,
        "v13_4_reference_step1000": 0.9371,
        "reproduce": "python scripts/v14_1_long_projection.py 5000",
    }, indent=2))
    print("RUN_DONE", json.dumps({
        "t_final": float(t_hist[-1]),
        "H_fine_final": float(H_f_hist[-1]),
        "orbits": phase_tot / (2.0 * np.pi),
    }))


if __name__ == "__main__":
    main()
