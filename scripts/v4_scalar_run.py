import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 50
print(f"N={N} L={L} h={h:.2f} dt={dt} steps={steps}")

# Gaussian scalar field (phi) and Pi = 0
s = minkowski_state(N)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z**2

for A in [1e-3, 1e-2]:
    print(f"\n== A = {A} ==")
    phi_f = A * np.exp(-R2 / 0.5)
    Pi_f  = np.zeros((N,N,N))
    s_run = s
    H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
        s_run.phi, s_run.gtilde, s_run.K, s_run.Atilde, h))))
    print(f"H0 = {H0:.4e}")
    t0 = time.time()
    for step in range(1, steps+1):
        s_run, phi_f, Pi_f = rk4_step_matter(s_run, phi_f, Pi_f, h, dt,
                                              V_func=lambda p: 0.5*p**2)
        if step % 10 == 0:
            H = float(np.mean(np.abs(hamiltonian_constraint_full(
                s_run.phi, s_run.gtilde, s_run.K, s_run.Atilde, h))))
            Kmax = np.max(np.abs(s_run.K))
            phimax = np.max(np.abs(s_run.phi))
            pi_max = np.max(np.abs(Pi_f))
            print(f"step={step:3d} H={H:.3e} Kmax={Kmax:.2e} "
                  f"phimax={phimax:.2e} Pimax={pi_max:.2e} "
                  f"finite={np.all(np.isfinite(s_run.phi)) and np.all(np.isfinite(Pi_f))} "
                  f"t={time.time()-t0:.0f}s")
