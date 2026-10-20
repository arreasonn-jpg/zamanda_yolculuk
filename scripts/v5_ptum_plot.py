import numpy as np, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_ptum import rk4_step_ptum
from zt005.matter.ptum_field import fuel_profile, m_eff_squared, psi_equilibrium
from zt005.matter.matter_ptum import ptum_sources

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0

Phi = fuel_profile(N, L, Phi0=3.0, sigma=4.0)
psi = psi_equilibrium(Phi, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Pi_Phi = np.zeros((N,N,N)); Pi_psi = np.zeros((N,N,N))

s = minkowski_state(N)
src0 = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                    V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                    lam_psi=lam_psi)

for step in range(steps):
    s, Phi, Pi_Phi, psi, Pi_psi = rk4_step_ptum(
        s, Phi, Pi_Phi, psi, Pi_psi, h, dt,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)

src1 = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                    V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                    lam_psi=lam_psi)

x = np.linspace(-L/2, L/2, N)
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# 1) rho along x-axis (y=z=0)
mid = N//2
axes[0,0].plot(x, src0["rho"][:, mid, mid], 'b-', label='t=0')
axes[0,0].plot(x, src1["rho"][:, mid, mid], 'r--', label=f't={dt*steps}')
axes[0,0].axhline(0, color='k', lw=0.5)
axes[0,0].set_xlabel('x'); axes[0,0].set_ylabel('rho')
axes[0,0].set_title('Energy density along x-axis')
axes[0,0].legend(); axes[0,0].grid(True)

# 2) rho 2D slice
im = axes[0,1].imshow(src1["rho"][:, :, mid], cmap='RdBu_r',
                       extent=[-L/2, L/2, -L/2, L/2], origin='lower')
plt.colorbar(im, ax=axes[0,1])
axes[0,1].set_title(f'rho 2D slice (t={dt*steps:.3f})')

# 3) negative cells count
neg_mask = src1["rho"] < 0
axes[0,2].imshow(neg_mask[:, :, mid].astype(int), cmap='Reds',
                  extent=[-L/2, L/2, -L/2, L/2], origin='lower')
axes[0,2].set_title(f'NEGATIVE cells: {np.sum(neg_mask)}/729')

# 4) K along x-axis
axes[1,0].plot(x, s.K[:, mid, mid], 'g-')
axes[1,0].axhline(0, color='k', lw=0.5)
axes[1,0].set_xlabel('x'); axes[1,0].set_ylabel('K')
axes[1,0].set_title('Extrinsic curvature K')
axes[1,0].grid(True)

# 5) phi metric perturbation
axes[1,1].plot(x, s.phi[:, mid, mid], 'm-')
axes[1,1].set_xlabel('x'); axes[1,1].set_ylabel('phi_conformal')
axes[1,1].set_title('Metric conformal factor')
axes[1,1].grid(True)

# 6) alpha (lapse)
axes[1,2].plot(x, s.alpha[:, mid, mid], 'c-')
axes[1,2].set_xlabel('x'); axes[1,2].set_ylabel('alpha')
axes[1,2].set_title('Lapse')
axes[1,2].grid(True)

plt.tight_layout()
plt.savefig("results/validated/v5_ptum_diagnostics.png", dpi=100)
print("Plot saved: results/validated/v5_ptum_diagnostics.png")
print(f"rho_min = {src1['rho'].min():.3e}")
print(f"neg_cells = {np.sum(neg_mask)}/729")
