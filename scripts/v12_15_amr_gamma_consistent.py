"""v12.15: constraint-consistent injection (Gamma rederivation).

Hypothesis: the fine-grid constraint violation injected by the coarse
BC (the long-standing v12.9-v12.13 limit) is reduced if the stored
fine Gamma^i in the injected band is rederived from the prolonged
gtilde with fine spacing (matching the RHS's own derivative operator)
instead of being prolongated from coarse stencils.

Setup identical to scripts/v12_9_amr_test.py (N_c=N_f=17, L_c=20,
L_f=10, N_sub=9, factor=2, dt=2e-3, eta=1.0, n_layers=4). Two runs:

  baseline     : injection as v12.9 (prolongated Gamma)
  gamma_cons   : Gamma in injected band rederived from fresh gtilde

Success criterion: H_f/H0_fine at step N < baseline's 1.0166
(v12.9/v12.11 trilinear reference at 1000 steps).
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
    prolongate_scalar, prolongate_tensor, prolongate_vector, inject_boundary)
from zt005.amr.bssn_amr_consistent import rederive_Gamma

N_c, L_c = 17, 20.0
N_f, L_f = 17, 10.0
N_sub = 9
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
DT = 2e-3


def extract_and_prolongate(s_c):
    lo = (N_c - N_sub) // 2
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


def _band_inject(fine, coarse, n_layers):
    out = fine.copy()
    for k in range(n_layers):
        out[k, ...] = coarse[k, ...]
        out[-1 - k, ...] = coarse[-1 - k, ...]
        out[:, k, ...] = coarse[:, k, ...]
        out[:, -1 - k, ...] = coarse[:, -1 - k, ...]
        out[:, :, k, ...] = coarse[:, :, k, ...]
        out[:, :, -1 - k, ...] = coarse[:, :, -1 - k, ...]
    return out


def state_inject(s_old, s_fresh, n_layers, gamma_consistent):
    fields = {
        "phi": _band_inject(s_old.phi, s_fresh.phi, n_layers),
        "gtilde": _band_inject(s_old.gtilde, s_fresh.gtilde, n_layers),
        "K": _band_inject(s_old.K, s_fresh.K, n_layers),
        "Atilde": _band_inject(s_old.Atilde, s_fresh.Atilde, n_layers),
        "alpha": _band_inject(s_old.alpha, s_fresh.alpha, n_layers),
        "beta": _band_inject(s_old.beta, s_fresh.beta, n_layers),
        "B": _band_inject(s_old.B, s_fresh.B, n_layers),
    }
    if gamma_consistent:
        Gam = np.array(s_old.Gamma, copy=True)
        Gam_f = rederive_Gamma(s_fresh.gtilde, h_f)
        for k in range(n_layers):
            Gam[k, :, :] = Gam_f[k, :, :]
            Gam[N_f - 1 - k, :, :] = Gam_f[N_f - 1 - k, :, :]
            Gam[:, k, :] = Gam_f[:, k, :]
            Gam[:, N_f - 1 - k, :] = Gam_f[:, N_f - 1 - k, :]
            Gam[:, :, k] = Gam_f[:, :, k]
            Gam[:, :, N_f - 1 - k] = Gam_f[:, :, N_f - 1 - k]
    else:
        Gam = _band_inject(s_old.Gamma, s_fresh.Gamma, n_layers)
    fields["Gamma"] = Gam
    return BSSNStateMP(**fields)


def run(mode, steps, log_every=200):
    s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
    s_f = extract_and_prolongate(s_c)
    gc = (mode == "gamma_cons")

    xc = np.linspace(-L_c / 2, L_c / 2, N_c)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    Rc = np.sqrt(Xc ** 2 + Yc ** 2 + Zc ** 2)
    lo, hi = (N_c - N_sub) // 2, (N_c - N_sub) // 2 + N_sub
    mask_c_sub = (Rc > 1.5)[lo:hi, lo:hi, lo:hi]
    xf = np.linspace(-L_f / 2, L_f / 2, N_f)
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    mask_f = np.sqrt(Xf ** 2 + Yf ** 2 + Zf ** 2) > 1.5

    H0_c_sub = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
        s_c.K[lo:hi, lo:hi, lo:hi],
        s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
    H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

    t_hist, H_c_hist, H_f_hist = [], [], []
    t0 = time.time()
    for step in range(1, steps + 1):
        s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        s_fresh = extract_and_prolongate(s_c)
        s_f = state_inject(s_f, s_fresh, 4, gc)
        s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=DT, width_cells=3,
                                  kappa=0.0, eta=1.0)
        if step % log_every == 0 or step == steps:
            H_c = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_c.phi[lo:hi, lo:hi, lo:hi],
                s_c.gtilde[lo:hi, lo:hi, lo:hi],
                s_c.K[lo:hi, lo:hi, lo:hi],
                s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
            H_f = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
            t_hist.append(step * DT)
            H_c_hist.append(H_c / H0_c_sub)
            H_f_hist.append(H_f / H0_f)
            fin = np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
            print(f"[{mode:>11s}] step={step:4d} "
                  f"H_c/H0={H_c / H0_c_sub:.4f} H_f/H0={H_f / H0_f:.4f} "
                  f"fin={fin} t={time.time() - t0:.0f}s", flush=True)
            if not fin:
                print("NON-FINITE"); break
    return {"mode": mode, "steps": steps,
            "H0_coarse_sub": H0_c_sub, "H0_fine": H0_f,
            "t": t_hist, "H_coarse_sub": H_c_hist, "H_fine": H_f_hist,
            "wall_s": time.time() - t0}


if __name__ == "__main__":
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    only = sys.argv[2] if len(sys.argv) > 2 else None
    modes = (only,) if only in ("baseline", "gamma_cons") else ("baseline", "gamma_cons")
    out = {"stage": "v12.15", "H_fine_reference_v12_9_step1000": 1.0166,
           "runs": []}
    for mode in modes:
        out["runs"].append(run(mode, steps))
    h_base = next((r["H_fine"][-1] for r in out["runs"]
                   if r["mode"] == "baseline"), None)
    gc_fin = next((r["H_fine"][-1] for r in out["runs"]
                   if r["mode"] == "gamma_cons"), None)
    out["H_baseline"] = h_base
    out["H_gamma_cons"] = gc_fin
    if h_base is not None and gc_fin is not None:
        if gc_fin < h_base:
            out["verdict"] = "gamma_consistent-BETTER"
        elif abs(gc_fin - h_base) < 0.002:
            out["verdict"] = "gamma_consistent-NEUTRAL"
        else:
            out["verdict"] = "gamma_consistent-WORSE"
    (base := Path("results/validated")).mkdir(exist_ok=True)
    for r in out["runs"]:
        np.savez(base / f"v12_15_amr_{r['mode']}.npz",
                 t=np.array(r["t"]), H_coarse_sub=np.array(r["H_coarse_sub"]),
                 H_fine=np.array(r["H_fine"]))
    (base / "v12_15_amr_gamma_consistent.json").write_text(json.dumps(out, indent=2))
    keys = ("H_baseline", "H_gamma_cons", "verdict")
    print(json.dumps({k: out[k] for k in keys if k in out}, indent=2))
