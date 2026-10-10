import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.sc_condensate import sc_condensate_profile, sc_condensate_potential

N = 9; L = 20.0; h = L/(N-1); dt = 1e-6; steps = 50
print(f"N={N} L={L} h={h:.2f} dt={dt} steps={steps}")

# Normalize: use small amplitude to keep numerics in linear regime
phi_f, Pi_f, sigma, omega, A_amp = sc_condensate_profile(
    N=N, L=L, M=1.0, rho_core=1e19, m_pair_MeV=10.0)
print(f"sigma = {sigma:.3e} m")
print(f"omega = {omega:.3e} 1/s")
print(f"A_amp = {A_amp:.3e}")

# Rescale phi to O(1e-3) to keep linear regime
scale = 1e-3 / A_amp
phi_f = phi_f * scale
print(f"phi_max after scale = {np.max(np.abs(phi_f)):.3e}")
print(f"phi_max original = {np.max(np.abs(phi_f))/scale:.3e}")

# Rescaled potential (small numbers)
V_func = sc_condensate_potential(m_pair_MeV=10.0, g_int=0.1)

s = minkowski_state(N)
t0 = time.time()
for step in range(1, steps+1):
    s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt, V_func=V_func)
    if step % 10 == 0:
        Kmax = np.max(np.abs(s.K))
        phimax = np.max(np.abs(phi_f))
        pimax = np.max(np.abs(Pi_f))
        fin = np.all(np.isfinite(s.phi)) and np.all(np.isfinite(phi_f)) and np.all(np.isfinite(Pi_f))
        print(f"step={step:3d} Kmax={Kmax:.3e} phimax={phimax:.3e} "
              f"Pimax={pimax:.3e} finite={fin} t={time.time()-t0:.0f}s")
