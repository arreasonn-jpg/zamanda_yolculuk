"""v12.11: AMR with cubic prolongation vs trilinear."""
import numpy as np, time, sys
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector,
    prolongate_scalar_cubic, prolongate_tensor_cubic, prolongate_vector_cubic,
)

N_c, L_c = 17, 20.0
N_sub = 9
L_f = 10.0
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 500
tag = sys.argv[2] if len(sys.argv) > 2 else "lin"
use_cubic = (sys.argv[3] if len(sys.argv) > 3 else "0") == "1"
factor = int(sys.argv[4]) if len(sys.argv) > 4 else 2

N_f = (N_sub - 1) * factor + 1
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
dt = 2e-3
order = "cubic" if use_cubic else "linear"
print(f"Coarse N={N_c} h={h_c:.4f} | Fine N={N_f} h={h_f:.4f} | {order} factor={factor}")

if use_cubic:
    ps, pt, pv = prolongate_scalar_cubic, prolongate_tensor_cubic, prolongate_vector_cubic
else:
    ps, pt, pv = prolongate_scalar, prolongate_tensor, prolongate_vector

def extract(s_c, N_c, N_sub, factor):
    lo = (N_c - N_sub) // 2
    hi = lo + N_sub
    return BSSNStateMP(
        phi=ps(s_c.phi[lo:hi, lo:hi, lo:hi], factor),
        gtilde=pt(s_c.gtilde[lo:hi, lo:hi, lo:hi], factor),
        K=ps(s_c.K[lo:hi, lo:hi, lo:hi], factor),
        Atilde=pt(s_c.Atilde[lo:hi, lo:hi, lo:hi], factor),
        Gamma=pv(s_c.Gamma[lo:hi, lo:hi, lo:hi], factor),
        alpha=ps(s_c.alpha[lo:hi, lo:hi, lo:hi], factor),
        beta=pv(s_c.beta[lo:hi, lo:hi, lo:hi], factor),
        B=pv(s_c.B[lo:hi, lo:hi, lo:hi], factor),
    )

def _inj(fine, coarse, n):
    out = fine.copy()
    for k in range(n):
        out[k, ...] = coarse[k, ...]; out[-1-k, ...] = coarse[-1-k, ...]
        out[:, k, ...] = coarse[:, k, ...]; out[:, -1-k, ...] = coarse[:, -1-k, ...]
        out[:, :, k, ...] = coarse[:, :, k, ...]; out[:, :, -1-k, ...] = coarse[:, :, -1-k, ...]
    return out

def inject(s_old, s_fresh, n=4):
    return BSSNStateMP(
        phi=_inj(s_old.phi, s_fresh.phi, n),
        gtilde=_inj(s_old.gtilde, s_fresh.gtilde, n),
        K=_inj(s_old.K, s_fresh.K, n),
        Atilde=_inj(s_old.Atilde, s_fresh.Atilde, n),
        Gamma=_inj(s_old.Gamma, s_fresh.Gamma, n),
        alpha=_inj(s_old.alpha, s_fresh.alpha, n),
        beta=_inj(s_old.beta, s_fresh.beta, n),
        B=_inj(s_old.B, s_fresh.B, n),
    )

s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
s_f = extract(s_c, N_c, N_sub, factor)

lo, hi = (N_c - N_sub) // 2, (N_c - N_sub) // 2 + N_sub
xc = np.linspace(-L_c/2, L_c/2, N_c)
Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
mask_c_sub = (np.sqrt(Xc**2 + Yc**2 + Zc**2) > 1.5)[lo:hi, lo:hi, lo:hi]
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
mask_f = np.sqrt(Xf**2 + Yf**2 + Zf**2) > 1.5

H0_c = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
    s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

t_hist, Hc_hist, Hf_hist = [], [], []
t0 = time.time()
for step in range(1, steps + 1):
    s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=dt, width_cells=3, kappa=0.0, eta=1.0)
    s_f_fresh = extract(s_c, N_c, N_sub, factor)
    s_f = inject(s_f, s_f_fresh, n=4)
    s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=dt, width_cells=3, kappa=0.0, eta=1.0)
    if step % 50 == 0:
        Hc = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
            s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
        Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
        t_hist.append(step*dt); Hc_hist.append(Hc/H0_c); Hf_hist.append(Hf/H0_f)
        fin = np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
        print(f"step={step:4d} H_c/H0={Hc/H0_c:.4f} H_f/H0={Hf/H0_f:.4f} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

np.savez(f"results/validated/v12_11_amr_{tag}.npz",
         t=np.array(t_hist), H_coarse=np.array(Hc_hist), H_fine=np.array(Hf_hist))
print(f"Saved v12_11_amr_{tag}")
