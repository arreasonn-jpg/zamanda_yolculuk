"""v12.12: 2x AMR long run + BH tracking (2000 steps)."""
import numpy as np, time, sys
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr_interp import prolongate_scalar, prolongate_tensor, prolongate_vector
from zt005.amr_bh_track import find_bh_centers, pair_separation

N_c, L_c, N_sub, L_f = 17, 20.0, 9, 10.0
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
tag = sys.argv[2] if len(sys.argv) > 2 else "long"
factor = 2
N_f = (N_sub - 1) * factor + 1
h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
dt = 2e-3

def extract(s_c):
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

def _inj(fine, coarse, n):
    out = fine.copy()
    for k in range(n):
        out[k, ...] = coarse[k, ...]; out[-1-k, ...] = coarse[-1-k, ...]
        out[:, k, ...] = coarse[:, k, ...]; out[:, -1-k, ...] = coarse[:, -1-k, ...]
        out[:, :, k, ...] = coarse[:, :, k, ...]; out[:, :, -1-k, ...] = coarse[:, :, -1-k, ...]
    return out

def inject(s_old, s_fresh, n=4):
    return BSSNStateMP(
        phi=_inj(s_old.phi, s_fresh.phi, n), gtilde=_inj(s_old.gtilde, s_fresh.gtilde, n),
        K=_inj(s_old.K, s_fresh.K, n), Atilde=_inj(s_old.Atilde, s_fresh.Atilde, n),
        Gamma=_inj(s_old.Gamma, s_fresh.Gamma, n), alpha=_inj(s_old.alpha, s_fresh.alpha, n),
        beta=_inj(s_old.beta, s_fresh.beta, n), B=_inj(s_old.B, s_fresh.B, n),
    )

s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
s_f = extract(s_c)
lo, hi = (N_c - N_sub) // 2, (N_c - N_sub) // 2 + N_sub
xc = np.linspace(-L_c/2, L_c/2, N_c)
Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
mask_c = (np.sqrt(Xc**2 + Yc**2 + Zc**2) > 1.5)[lo:hi, lo:hi, lo:hi]
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
mask_f = np.sqrt(Xf**2 + Yf**2 + Zf**2) > 1.5
H0_c = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
    s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c])))
H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

t_hist, Hc_hist, Hf_hist, xl_hist, xr_hist, sep_hist = [], [], [], [], [], []
t0 = time.time()
for step in range(1, steps + 1):
    s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=dt, width_cells=3, kappa=0.0, eta=1.0)
    s_f_fresh = extract(s_c)
    s_f = inject(s_f, s_f_fresh, n=4)
    s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=dt, width_cells=3, kappa=0.0, eta=1.0)
    if step % 100 == 0:
        xl, xr, _ = find_bh_centers(s_c, L_c, N_c)
        sep = pair_separation(xl, xr)
        Hc = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
            s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c])))
        Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
        t_hist.append(step*dt); Hc_hist.append(Hc/H0_c); Hf_hist.append(Hf/H0_f)
        xl_hist.append(xl); xr_hist.append(xr); sep_hist.append(sep)
        fin = np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
        print(f"step={step:5d} H_c/H0={Hc/H0_c:.4f} H_f/H0={Hf/H0_f:.4f} "
              f"sep={sep:.2f} xl={xl:.2f} xr={xr:.2f} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

np.savez(f"results/validated/v12_12_amr_{tag}.npz",
         t=np.array(t_hist), H_coarse=np.array(Hc_hist), H_fine=np.array(Hf_hist),
         xl=np.array(xl_hist), xr=np.array(xr_hist), sep=np.array(sep_hist))
print(f"Saved v12_12_amr_{tag}")
