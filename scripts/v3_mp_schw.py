import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_puncture import schwarzschild_puncture_state
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 40.0; M = 1.0; dt = 5e-4; steps = 50
h = L/(N-1)
print(f"N={N} L={L} h={h:.2f} dt={dt} steps={steps}")

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
    if step % 10 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        n = float(np.mean(np.abs(H[mask])))
        norms.append(n)
        print(f"step={step:3d} H_ext={n:.3e} growth={n/H0:.4f} "
              f"alpha_min={s.alpha.min():.4f} Bmax={np.max(np.abs(s.B)):.2e} "
              f"t={time.time()-t0:.0f}s")
print(f"final growth = {norms[-1]/norms[0]:.4f}")
