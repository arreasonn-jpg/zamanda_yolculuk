"""v12.10: AMR with configurable refinement factor and sub-cycling.

Usage: v12_10_amr_flex.py STEPS TAG FACTOR SUBCYCLE
  FACTOR=2,4,... refinement (h_fine = h_coarse / FACTOR)
  SUBCYCLE=1,2,... fine steps per coarse step (dt_fine = dt_coarse / SUBCYCLE)
"""
import numpy as np, time, sys
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector, inject_boundary,
)

N_c, L_c = 17, 20.0
N_sub = 9
L_f = 10.0

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 500
tag = sys.argv[2] if len(sys.argv) > 2 else "v1"
factor = int(sys.argv[3]) if len(sys.argv) > 3 else 2
subcycle = int(sys.argv[4]) if len(sys.argv) > 4 else 1

N_f = (N_sub - 1) * factor + 1
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
dt_c = 2e-3
dt_f = dt_c / subcycle

print(f"Coarse N={N_c} h={h_c:.4f} dt={dt_c}")
print(f"Fine N={N_f} h={h_f:.4f} dt={dt_f} refinement={h_c/h_f:.2f}x subcycle={subcycle}")

def extract_and_prolongate(s_c, N_c, N_sub, factor=2):
    lo = (N_c - N_sub) // 2
    hi = lo + N_sub
    return BSSNStateMP(
        phi=prolongate_scalar(s_c.phi[lo:hi, lo:hi, lo:hi], factor),
        gtilde=prolongate_tensor(s_c.gtilde[lo:hi, lo:hi, lo:hi], factor),
        K=prolongate_scalar(s_c.K[lo:hi, lo:hi, lo:hi], factor),
        Atilde=prolongate_tensor(s_c.Atilde[lo:hi, lo:hi, lo:hi], factor),
        Gamma=prolongate_vector(s_c.Gamma[lo:hi, lo:hi, lo:hi], factor),
        alpha=prolongate_scalar(s_c.alpha[lo:hi, lo:hi, lo:hi], factor),
        beta=prolongate_vector(s_c.beta[lo:hi, lo:hi, lo:hi], factor),
        B=prolongate_vector(s_c.B[lo:hi, lo:hi, lo:hi], factor),
    )

def _inject_nd(fine, coarse, n_layers):
    out = fine.copy()
    for k in range(n_layers):
        out[k, ...] = coarse[k, ...]
        out[-1-k, ...] = coarse[-1-k, ...]
        out[:, k, ...] = coarse[:, k, ...]
        out[:, -1-k, ...] = coarse[:, -1-k, ...]
        out[:, :, k, ...] = coarse[:, :, k, ...]
        out[:, :, -1-k, ...] = coarse[:, :, -1-k, ...]
    return out

def inject_full(s_old, s_fresh, n_layers=4):
    return BSSNStateMP(
        phi=_inject_nd(s_old.phi, s_fresh.phi, n_layers),
        gtilde=_inject_nd(s_old.gtilde, s_fresh.gtilde, n_layers),
        K=_inject_nd(s_old.K, s_fresh.K, n_layers),
        Atilde=_inject_nd(s_old.Atilde, s_fresh.Atilde, n_layers),
        Gamma=_inject_nd(s_old.Gamma, s_fresh.Gamma, n_layers),
        alpha=_inject_nd(s_old.alpha, s_fresh.alpha, n_layers),
        beta=_inject_nd(s_old.beta, s_fresh.beta, n_layers),
        B=_inject_nd(s_old.B, s_fresh.B, n_layers),
    )

s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
s_f = extract_and_prolongate(s_c, N_c, N_sub, factor=factor)

lo = (N_c - N_sub) // 2
hi = lo + N_sub
xc = np.linspace(-L_c/2, L_c/2, N_c)
Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
Rc = np.sqrt(Xc**2 + Yc**2 + Zc**2)
mask_c_sub = (Rc > 1.5)[lo:hi, lo:hi, lo:hi]

xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
mask_f = np.sqrt(Xf**2 + Yf**2 + Zf**2) > 1.5

H0_c_sub = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
    s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
print(f"H0_c_sub={H0_c_sub:.4e} H0_f={H0_f:.4e}")

t_hist, H_c_hist, H_f_hist = [], [], []
t0 = time.time()

for step in range(1, steps + 1):
    s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=dt_c, width_cells=3, kappa=0.0, eta=1.0)
    s_f_fresh = extract_and_prolongate(s_c, N_c, N_sub, factor=factor)
    for sub in range(subcycle):
        s_f = inject_full(s_f, s_f_fresh, n_layers=4)
        s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=dt_f, width_cells=3, kappa=0.0, eta=1.0)
    if step % 50 == 0:
        H_c_sub = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
            s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
        H_f = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
        t_hist.append(step * dt_c)
        H_c_hist.append(H_c_sub / H0_c_sub)
        H_f_hist.append(H_f / H0_f)
        fin = np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
        print(f"step={step:4d} H_c/H0={H_c_sub/H0_c_sub:.4f} H_f/H0={H_f/H0_f:.4f} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

np.savez(f"results/validated/v12_10_amr_{tag}.npz",
         t=np.array(t_hist), H_coarse=np.array(H_c_hist), H_fine=np.array(H_f_hist))
print(f"Saved v12_10_amr_{tag}")
