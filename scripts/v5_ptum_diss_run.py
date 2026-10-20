import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum_diss import rk4_step_ptum_diss
from zt005.matter.ptum_field import fuel_profile, m_eff_squared, psi_equilibrium
from zt005.matter.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 2e-4; steps = 200
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
eps_KO = 0.5
print(f"N={N} dt={dt} steps={steps} eps_KO={eps_KO}")

Phi = fuel_profile(N, L, Phi0=1.5, sigma=4.0)
psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Pi_Phi = np.zeros((N,N,N)); Pi_psi = np.zeros((N,N,N))

s = minkowski_state(N)
t0 = time.time()
for step in range(1, steps+1):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum_diss(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi,
        eps_KO=eps_KO)
    if step % 50 == 0:
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                           V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
        n_neg = int(np.sum(src["rho"] < 0))
        Kmax = float(np.max(np.abs(s.K)))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
               and np.all(np.isfinite(psi)))
        print(f"step={step:4d} rho_min={src['rho'].min():.3e} "
              f"neg={n_neg:3d}/729 Kmax={Kmax:.3e} finite={fin} "
              f"t={time.time()-t0:.0f}s")
