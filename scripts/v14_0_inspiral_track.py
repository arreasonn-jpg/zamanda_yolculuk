"""v14.0: instrumented long inspiral run (moving-patch AMR, no projection).

The first stage of the full-inspiral marathon. The v12.9-validated
fixed-patch 2x AMR harness (N_c = 17, window 9, factor 2, dt = 2e-3)
runs with the v13.0 COM-driven recentering ENABLED and full
orbital instrumentation added:

- per-side BH centers (x, y) via amr_bh_track.peak_centroids on the
  coarse grid (2D weighted centroids; a pure x-axis pair would be
  invisible for an orbit);
- separation d(t), orbital angle phi(t) (unwrapped), angular
  frequency omega(t) at every log point;
- Newtonian reference omega_N = sqrt(M / d^3);
- the standard H_c/H0 and H_f/H0 curves and the v13.0 recenter log.

The Bowen-York binary already carries orbital tangential momenta
(P = m v_orb, v_orb = sqrt(M/(4d))); previous 1000-step runs only
covered t = 2 (~4% of an orbit at d = 4), so the orbit was never
resolved. This run targets 10000 steps (t = 20, ~0.4 orbit).

No projection is applied: this run measures the baseline physics of
the inspiral, on top of the already-recorded AMR defect chain.
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
from zt005.amr.amr_moving import target_lo as target_lo_impl
from zt005.amr.amr_moving import fine_grid_mask
from zt005.amr.amr_bh_track import peak_centroids, alpha_pair

N_c, L_c = 17, 20.0
N_f, L_f = 17, 10.0
N_sub = 9
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
DT = 2e-3
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
M_TOTAL = 1.0
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
    s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5,
                                                  4.0))
    lo = target_lo_impl(0.0, h_c, N_c, N_sub)   # initial: centered
    s_f = extract_window(s_c, lo)
    mask_f = fine_mask()

    hi_c = CENTER_IDX + N_sub
    xc = np.linspace(-L_c / 2, L_c / 2, N_c)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    mask_c_sub = (np.sqrt(Xc ** 2 + Yc ** 2 + Zc ** 2)
                  > 1.5)[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                         CENTER_IDX:hi_c]
    H0_c = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_c.phi[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        s_c.gtilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                   CENTER_IDX:hi_c],
        s_c.K[CENTER_IDX:hi_c, CENTER_IDX:hi_c, CENTER_IDX:hi_c],
        s_c.Atilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                   CENTER_IDX:hi_c], h_c)[mask_c_sub])))
    H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

    t_hist, H_c_hist, H_f_hist = [], [], []
    d_hist, phase_hist, omega_hist = [], [], []
    centers_hist = []
    # alpha-funnel meter series (second, independent locator; see
    # amr_bh_track.alpha_pair for the calibration lesson)
    d_al_hist = []
    th_al_hist = []
    amin_hist = []
    recenter_steps = []
    phase_prev = None
    phase_tot = 0.0
    d0_ref = None
    t0 = time.time()
    for step in range(1, STEPS + 1):
        s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        # BH meter on the FINE window (h = 0.625, 2x the coarse
        # resolution; both punctures stay inside while the equal-mass
        # COM holds the window centered). Window-frame coordinates
        # shift back to global by off = x_c[lo] + L_f/2.
        # Two independent locators are recorded:
        #  - peak_centroids: the |Atilde|^2 bulge meter (calibration
        #    lesson: BY bulges are momentum rings at r > 0 around the
        #    punctures, so positions are ring-biased);
        #  - alpha_pair: the lapse-basin meter (calibration lesson:
        #    no funnel on the coarse grid within 6 M; observable at
        #    fine resolution only if/once it forms).
        (cl_a, cr_a, _) = peak_centroids(s_f, L_f)
        (cl_g, cr_g, depths_g) = alpha_pair(s_f, L_f)
        off = xc[lo] + L_f / 2.0
        xl, yl = cl_a[0] + off, cl_a[1] + off
        xr, yr = cr_a[0] + off, cr_a[1] + off
        lo_new = target_lo_impl(0.5 * (xl + xr), h_c, N_c, N_sub)
        if lo_new != lo:
            lo = lo_new
            s_f = extract_window(s_c, lo)
            recenter_steps.append(step)
            mask_f = fine_mask()
        s_fresh = extract_window(s_c, lo)
        s_f = inject_full(s_f, s_fresh, n_layers=4)
        s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        if step % 50 == 0 or step == STEPS:
            t = step * DT
            d = float(np.hypot(xr - xl, yr - yl))
            d_al = float(np.hypot(cr_g[0] - cl_g[0],
                                  cr_g[1] - cl_g[1]))
            th_al = float(np.arctan2(cr_g[1] - cl_g[1],
                                     cr_g[0] - cl_g[0]))
            amin_hist.append((float(depths_g[0]),
                              float(depths_g[1])))
            d_al_hist.append(d_al)
            th_al_hist.append(th_al)
            if d0_ref is None:
                d0_ref = d
            ang = float(np.arctan2(yr - yl, xr - xl))
            if phase_prev is None:
                dphase = 0.0
            else:
                dphase = ang - phase_prev
                while dphase > np.pi:
                    dphase -= 2.0 * np.pi
                while dphase < -np.pi:
                    dphase += 2.0 * np.pi
                phase_tot += dphase
                omega_hist.append(dphase / (t - t_hist[-1]))
            phase_prev = ang
            t_hist.append(t)
            phase_hist.append(phase_tot)
            d_hist.append(d)
            centers_hist.append(((xl, yl), (xr, yr)))
            om_last = omega_hist[-1] if omega_hist else float("nan")
            H_c = float(np.mean(np.abs(
                hamiltonian_constraint_full(
                    s_c.phi[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                            CENTER_IDX:hi_c],
                    s_c.gtilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                               CENTER_IDX:hi_c],
                    s_c.K[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                          CENTER_IDX:hi_c],
                    s_c.Atilde[CENTER_IDX:hi_c, CENTER_IDX:hi_c,
                               CENTER_IDX:hi_c],
                    h_c)[mask_c_sub])))
            H_f = float(np.mean(np.abs(
                hamiltonian_constraint_full(
                    s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde,
                    h_f)[mask_f])))
            H_c_hist.append(H_c / H0_c)
            H_f_hist.append(H_f / H0_f)
            fin = np.all(np.isfinite(s_c.phi)) \
                and np.all(np.isfinite(s_f.phi))
            print(f"[v14 insp] t={t:6.3f} d={d:.4f} "
                  f"phase={phase_tot:+.4f} "
                  f"omega={om_last:.4f} "
                  f"d_al={d_al:.4f} "
                  f"a_min={min(depths_g):.3f} "
                  f"H_c/H0={H_c / H0_c:.4f} H_f/H0={H_f / H0_f:.4f} "
                  f"rec={len(recenter_steps)} fin={fin} "
                  f"t WALL={time.time() - t0:.0f}s", flush=True)
            if not fin:
                print("NON-FINITE")
                break
            if step % 1000 == 0:
                base_ckpt = Path("results/validated")
                np.savez(base_ckpt / "v14_0_inspiral_track.npz",
                         t=np.array(t_hist),
                         H_coarse_sub=np.array(H_c_hist),
                         H_fine=np.array(H_f_hist),
                         d=np.array(d_hist),
                         phase=np.array(phase_hist),
                         omega=np.array(omega_hist),
                         d_alpha=np.array(d_al_hist),
                         theta_alpha=np.array(th_al_hist),
                         alpha_min=np.array(amin_hist))
                print(f"[checkpoint] v14.0 step={step}", flush=True)

    base = Path("results/validated")
    cl_arr = np.array([c[0] for c in centers_hist])   # (n, 2)
    cr_arr = np.array([c[1] for c in centers_hist])   # (n, 2)
    np.savez(base / "v14_0_inspiral_track.npz",
             t=np.array(t_hist), H_coarse_sub=np.array(H_c_hist),
             H_fine=np.array(H_f_hist), d=np.array(d_hist),
             phase=np.array(phase_hist),
             omega=np.array(omega_hist),
             d_alpha=np.array(d_al_hist),
             theta_alpha=np.array(th_al_hist),
             alpha_min=np.array(amin_hist),
             cl_x=cl_arr[:, 0], cl_y=cl_arr[:, 1],
             cr_x=cr_arr[:, 0], cr_y=cr_arr[:, 1])
    omega_N_ref = float(np.sqrt(M_TOTAL / d0_ref ** 3))
    mean_om = float(np.mean(omega_hist)) if omega_hist else 0.0
    (base / "v14_0_inspiral_track.json").write_text(json.dumps({
        "stage": "v14.0",
        "steps": STEPS,
        "t_final": t_hist[-1],
        "d0": d0_ref,
        "d_final": d_hist[-1],
        "d_ratio": d_hist[-1] / d0_ref,
        "phase_total_rad": phase_tot,
        "orbits": phase_tot / (2.0 * np.pi),
        "omega_measured_mean": mean_om,
        "omega_newtonian_ref": omega_N_ref,
        "omega_ratio": mean_om / omega_N_ref,
        "d_alpha_final": d_al_hist[-1],
        "alpha_min_final": float(min(amin_hist[-1])),
        "recenter_steps": recenter_steps,
        "recenters": len(recenter_steps),
        "H_fine_final": H_f_hist[-1],
        "H_coarse_final": H_c_hist[-1],
        "baseline_AMR_reference_step1000": 1.0756,
        "reproduce": "python scripts/v14_0_inspiral_track.py 10000",
    }, indent=2))
    print("RUN_DONE", json.dumps({
        "d_ratio": d_hist[-1] / d0_ref,
        "orbits": phase_tot / (2.0 * np.pi),
        "omega_ratio": mean_om / omega_N_ref,
        "recenters": len(recenter_steps)}))


if __name__ == "__main__":
    main()
