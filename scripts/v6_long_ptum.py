import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.ptum_field import fuel_profile, psi_equilibrium
from zt005.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 2e-4
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
A = 3.0

s = minkowski_state(N)
Phi = fuel_profile(N, L, Phi0=A, sigma=4.0)
psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Z = np.zeros((N,N,N))
Pi_Phi = Z.copy(); Pi_psi = Z.copy()

n_neg_history = []
print(f"LONG PTUM RUN A={A}, dt={dt}, target 1000 steps")
t0 = time.time()
for step in range(1, 1001):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    if step % 100 == 0:
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                           V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
        rho_min = float(src["rho"].min())
        n_neg = int(np.sum(src["rho"] < 0))
        Phi_max = float(Phi.max())
        Kmax = float(np.max(np.abs(s.K)))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
               and np.all(np.isfinite(psi)))
        n_neg_history.append((step, rho_min, n_neg, Phi_max, Kmax))
        print(f"step={step:4d} rho_min={rho_min:.4e} neg={n_neg:3d} "
              f"Phimax={Phi_max:.3e} Kmax={Kmax:.3e} finite={fin} "
              f"t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE — stopping")
            break

print(f"\nNegative energy lifetime:")
first_neg = next((s for s,_,n,_,_ in n_neg_history if n > 0), None)
last_neg = next((s for s,_,n,_,_ in reversed(n_neg_history) if n > 0), None)
print(f"  first non-zero: step {first_neg}")
print(f"  last non-zero:  step {last_neg}")
