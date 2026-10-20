"""v12.9-alpha2: fixed-patch AMR test.

Two grids in lockstep:
- Coarse: N=17, L=20, h=1.25 (full domain)
- Fine:   N=17, L=10, h=0.625 (central sub-domain, prolonged from coarse)

Per step:
1. Coarse RK4 step (full domain)
2. Extract coarse central N_sub=9^3 sub-region
3. Prolongate factor=2 to fine N=17 grid
4. Inject prolonged coarse into fine boundary (n_layers=4)
5. Fine RK4 step
6. Track H on both, normalized to initial

Goal: H_fine/H0_fine << H_coarse_sub/H0_coarse_sub -> AMR helps.
"""
import numpy as np, time, sys
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.amr.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector, inject_boundary,
)

N_c, L_c = 17, 20.0
N_f, L_f = 17, 10.0
N_sub = 9
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
tag = sys.argv[2] if len(sys.argv) > 2 else "amr1"

h_c = L_c / (N_c - 1)
h_f = L_f / (N_f - 1)
print(f"Coarse N={N_c} L={L_c} h={h_c:.4f}")
print(f"Fine   N={N_f} L={L_f} h={h_f:.4f}  refinement={h_c/h_f:.2f}x")
print(f"Sub-extract {N_sub}^3 from coarse -> {N_f}^3 fine")

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
    """Generic N-d inject_boundary preserving trailing dims."""
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
s_f = extract_and_prolongate(s_c, N_c, N_sub, factor=2)

xc = np.linspace(-L_c/2, L_c/2, N_c)
Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
Rc = np.sqrt(Xc**2 + Yc**2 + Zc**2)
lo = (N_c - N_sub) // 2
hi = lo + N_sub
mask_c_full = Rc > 1.5
mask_c_sub = mask_c_full[lo:hi, lo:hi, lo:hi]

xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
Rf = np.sqrt(Xf**2 + Yf**2 + Zf**2)
mask_f = Rf > 1.5

H0_c_sub = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
    s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
H0_f = float(np.mean(np.abs(hamiltonian_constraint_full(
    s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))

print(f"H0_coarse_sub = {H0_c_sub:.4e}")
print(f"H0_fine       = {H0_f:.4e}")

t_hist, H_c_hist, H_f_hist = [], [], []
t0 = time.time()

for step in range(1, steps + 1):
    s_c = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=2e-3, width_cells=3, kappa=0.0, eta=1.0)
    s_f_fresh = extract_and_prolongate(s_c, N_c, N_sub, factor=2)
    s_f = inject_full(s_f, s_f_fresh, n_layers=4)
    s_f = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=2e-3, width_cells=3, kappa=0.0, eta=1.0)
    if step % 100 == 0:
        H_c_sub = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_c.phi[lo:hi, lo:hi, lo:hi], s_c.gtilde[lo:hi, lo:hi, lo:hi],
            s_c.K[lo:hi, lo:hi, lo:hi], s_c.Atilde[lo:hi, lo:hi, lo:hi], h_c)[mask_c_sub])))
        H_f = float(np.mean(np.abs(hamiltonian_constraint_full(
            s_f.phi, s_f.gtilde, s_f.K, s_f.Atilde, h_f)[mask_f])))
        t_hist.append(step * 2e-3)
        H_c_hist.append(H_c_sub / H0_c_sub)
        H_f_hist.append(H_f / H0_f)
        fin = np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
        print(f"step={step:4d} H_c_sub/H0={H_c_sub/H0_c_sub:.4f} "
              f"H_f/H0={H_f/H0_f:.4f} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

np.savez(f"results/validated/v12_9_amr_{tag}.npz",
         t=np.array(t_hist),
         H_coarse_sub=np.array(H_c_hist),
         H_fine=np.array(H_f_hist))
print(f"Saved v12_9_amr_{tag}")
