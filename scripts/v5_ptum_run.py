import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.matter.ptum_field import fuel_profile, Phi_critical, m_eff_squared
from zt005.matter.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 500
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
print(f"N={N} L={L} dt={dt} steps={steps}")

Phi_c = Phi_critical(alpha_p, Eth)
print(f"Phi_c = {Phi_c}")

# Fuel: activated in wide region
Phi = fuel_profile(N, L, Phi0=3.0, sigma=4.0)

# Seed psi with small Gaussian (mimics quantum fluctuation)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z**2
psi_seed = 1e-3 * np.exp(-R2 / 8.0)
psi = psi_seed.copy()
Pi_Phi = np.zeros((N,N,N))
Pi_psi = np.zeros((N,N,N))

meff = m_eff_squared(Phi, m2, alpha_p, Eth)
print(f"Activated cells: {np.sum(meff < 0)}/{Phi.size}")
print(f"psi_seed max = {psi.max():.3e}")

s = minkowski_state(N)
t0 = time.time()
for step in range(1, steps+1):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    if step % 50 == 0:
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                           V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                           lam_psi=lam_psi)
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
               and np.all(np.isfinite(psi)))
        print(f"step={step:3d} rho_min={src['rho'].min():.3e} "
              f"rho_max={src['rho'].max():.3e} psimax={np.max(np.abs(psi)):.3e} Pipsimax={np.max(np.abs(Pi_psi)):.3e} "
              f"Kmax={np.max(np.abs(s.K)):.3e} finite={fin} t={time.time()-t0:.0f}s")
