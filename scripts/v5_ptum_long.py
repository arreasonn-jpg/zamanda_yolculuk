import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.ptum_field import fuel_profile, m_eff_squared, psi_equilibrium
from zt005.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 2e-4; steps = 1000
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
print(f"N={N} dt={dt} steps={steps} total_t={dt*steps}")

Phi = fuel_profile(N, L, Phi0=3.0, sigma=4.0)
meff = m_eff_squared(Phi, m2, alpha_p, Eth)
psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Pi_Phi = np.zeros((N,N,N)); Pi_psi = np.zeros((N,N,N))

s = minkowski_state(N)
t0 = time.time()
norms = []
for step in range(1, steps+1):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    if step % 100 == 0:
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                           V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
        n_neg = int(np.sum(src["rho"] < 0))
        rho_min = float(src["rho"].min())
        Kmax = float(np.max(np.abs(s.K)))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
               and np.all(np.isfinite(psi)))
        norms.append((step, rho_min, n_neg, Kmax))
        print(f"step={step:4d} rho_min={rho_min:.4e} neg={n_neg:3d}/729 "
              f"Kmax={Kmax:.3e} finite={fin} t={time.time()-t0:.0f}s")

print(f"\nRho_min: initial={norms[0][1]:.4e} final={norms[-1][1]:.4e}")
print(f"Neg cells: initial={norms[0][2]} final={norms[-1][2]}")
print(f"Kmax: initial={norms[0][3]:.3e} final={norms[-1][3]:.3e}")
