import numpy as np, time
from zt005.sc_full_ic import sc_full_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; M = 1.0; rc = 2.0; dt = 5e-4; steps = 50
h = L/(N-1)
print(f"N={N} L={L} M={M} rc={rc} dt={dt} steps={steps}")

s, phi_f, Pi_f, V_func = sc_full_state(N=N, L=L, M=M, rc=rc,
                                        phi_amp_target=1e-4)
print(f"phi_max = {np.max(np.abs(phi_f)):.3e}")
print(f"alpha_min = {s.alpha.min():.4f}")

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
mask = R > 3.0*rc

H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
print(f"H0 = {H0:.4e}")

t0 = time.time()
for step in range(1, steps+1):
    s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt, V_func=V_func)
    if step % 10 == 0:
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        nH = float(np.mean(np.abs(H[mask])))
        phimax = np.max(np.abs(phi_f))
        pimax = np.max(np.abs(Pi_f))
        Kmax = np.max(np.abs(s.K))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(phi_f))
               and np.all(np.isfinite(Pi_f)))
        print(f"step={step:3d} H={nH:.3e} growth={nH/H0:.4f} "
              f"phimax={phimax:.2e} Pimax={pimax:.2e} Kmax={Kmax:.2e} "
              f"alpha_min={s.alpha.min():.4f} finite={fin} "
              f"t={time.time()-t0:.0f}s")
