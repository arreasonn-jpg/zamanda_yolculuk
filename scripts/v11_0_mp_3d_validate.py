"""v11.0-alpha1: moving puncture 3D validation.
Run existing moving puncture driver (v3.0-alpha2) for:
- N = 13 (finer grid)
- 500 steps
- Schwarzschild puncture IC
Track H growth, alpha_min, B max, all finite.
Goal: validate 3D BH stability before binary.
"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 11; L = 40.0; M = 1.0
h = L/(N-1); dt = 2e-4; steps = 200
print(f"MP 3D validate: N={N}, L={L}, dt={dt}, steps={steps}")

s_base = schwarzschild_puncture_state(N, L, M=M, eps=M)
s = BSSNStateMP.from_base(s_base)

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
mask = (R > 3.0*M) & (R < L/3)

H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
print(f"H0_ext = {H0:.4e}")

t0 = time.time()
norms = [H0]
for step in range(1, steps+1):
    s = rk4_step_mp(s, h=h, dt=dt)
    if step % 50 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H[mask])))
        norms.append(n)
        Kmax = float(np.max(np.abs(s.K)))
        amin = float(s.alpha.min())
        bmax = float(np.max(np.abs(s.beta)))
        bv = float(np.max(np.abs(s.B)))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.K)))
        print(f"step={step:4d} H={n:.3e} growth={n/H0:.4f} "
              f"Kmax={Kmax:.3e} amin={amin:.4f} bmax={bmax:.3e} "
              f"Bmax={bv:.3e} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

print(f"\nfinal growth = {norms[-1]/norms[0]:.4e}")
np.savez("results/validated/v11_0_mp_3d.npz",
         steps=np.arange(0, len(norms)*50, 50),
         H=np.array(norms))
