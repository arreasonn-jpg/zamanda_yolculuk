import numpy as np, time
from zt005.kerr_sc_full_ic import kerr_sc_full_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4
total_steps = 500
rc = 2.0; phi_amp_target = 1e-4
a = 0.9  # worst case (near-extremal)
print(f"KERR+SC LONG: N={N} a={a} dt={dt} steps={total_steps}")

s, phi_f, Pi_f, V_func = kerr_sc_full_state(
    N=N, L=L, M=1.0, a=a, rc=rc, phi_amp_target=phi_amp_target)

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
mask = (R > 3.0*rc) & (R < L/3)

H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
print(f"H0 = {H0:.4e}")

t0 = time.time()
history = []
for step in range(1, total_steps+1):
    s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt, V_func=V_func)
    if step % 50 == 0:
        H = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        phimax = float(np.max(np.abs(phi_f)))
        Kmax = float(np.max(np.abs(s.K)))
        amin = float(s.alpha.min())
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(phi_f))
               and np.all(np.isfinite(Pi_f)))
        history.append((step, H, phimax, Kmax, amin))
        print(f"step={step:4d} H={H:.4e} growth={H/H0:.4f} "
              f"phimax={phimax:.3e} Kmax={Kmax:.3e} alpha_min={amin:.4f} "
              f"fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE — stopping")
            break

np.savez("results/validated/v6_kerr_sc_500.npz",
         steps=np.array([h[0] for h in history]),
         H=np.array([h[1] for h in history]),
         Kmax=np.array([h[3] for h in history]))
print("Saved: results/validated/v6_kerr_sc_500.npz")
