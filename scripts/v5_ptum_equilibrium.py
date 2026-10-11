import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.ptum_field import (fuel_profile, Phi_critical, m_eff_squared,
                                psi_equilibrium)
from zt005.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
print(f"N={N} L={L} dt={dt} steps={steps}")

Phi = fuel_profile(N, L, Phi0=3.0, sigma=4.0)
meff = m_eff_squared(Phi, m2, alpha_p, Eth)
print(f"Activated cells: {np.sum(meff < 0)}/{Phi.size}")

# START AT EQUILIBRIUM (SC condensate already formed)
psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
print(f"psi_eq max = {psi.max():.4f}")

Pi_Phi = np.zeros((N,N,N))
Pi_psi = np.zeros((N,N,N))

s = minkowski_state(N)
src0 = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                    V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                    lam_psi=lam_psi)
print(f"Initial rho: min={src0['rho'].min():.3e} max={src0['rho'].max():.3e}")
print(f"NEGATIVE CELLS: {np.sum(src0['rho'] < 0)}/{src0['rho'].size}")

t0 = time.time()
for step in range(1, steps+1):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    if step % 20 == 0:
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                           V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
        n_neg = np.sum(src["rho"] < 0)
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
               and np.all(np.isfinite(psi)))
        print(f"step={step:3d} rho_min={src['rho'].min():.3e} "
              f"neg_cells={n_neg}/{src['rho'].size} "
              f"psimax={np.max(np.abs(psi)):.3e} Kmax={np.max(np.abs(s.K)):.3e} "
              f"finite={fin} t={time.time()-t0:.0f}s")
