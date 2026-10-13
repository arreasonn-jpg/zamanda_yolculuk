"""v8.1-alpha3: LQC bounce with signed H."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.lqc_scalar_evolution import evolve_lqc_scalar

RHO_C = 1.0
G_EFF = 3.0 / (8.0 * np.pi)

# Start on contracting branch with rho close to rho_c
V0 = 0.5
V = lambda phi: V0 * np.exp(-phi**2 / 4.0)
dV = lambda phi: -phi/2.0 * V0 * np.exp(-phi**2 / 4.0)

# Choose phi_dot such that rho0 ~ 0.85 * rho_c
phi0 = 0.0
rho_target = 0.85 * RHO_C
phi_dot0 = -np.sqrt(max(2*(rho_target - V0), 0.0))

print(f"Bounce: phi0={phi0}, phi_dot0={phi_dot0:.4f}, rho0~{rho_target}")

out = evolve_lqc_scalar(phi0=phi0, phi_dot0=phi_dot0, V_func=V, dV_func=dV,
                         rho_c=RHO_C, dt=1e-3, n_steps=30000, G_eff=G_EFF,
                         H0_sign=-1)

rho = out["rho"]; H = out["H"]; a = out["a"]; t = out["t"]; phi = out["phi"]
print(f"\nphi: {phi[0]:.4f} -> {phi[-1]:.4f}")
print(f"rho: {rho[0]:.4e} -> max {rho.max():.4e} -> end {rho[-1]:.4e}")
print(f"rho_max/rho_c = {rho.max()/RHO_C:.4f}")
# Detect H sign change (bounce)
sign_flip = np.where(np.diff(np.sign(H)))[0]
if len(sign_flip) > 0:
    t_bounce = t[sign_flip[0]+1]
    print(f"BOUNCE detected at t = {t_bounce:.4f}")
else:
    print("No bounce detected (H stays same sign)")
print(f"H: init {H[0]:.4e}, end {H[-1]:.4e}")
print(f"a: init {a[0]:.4e}, min {a.min():.4e}, end {a[-1]:.4e}")

np.savez("results/validated/v8_1_bounce.npz",
         t=t, phi=phi, rho=rho, H=H, a=a)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes[0,0].plot(t, phi, 'b-'); axes[0,0].set_xlabel('t')
axes[0,0].set_ylabel('phi'); axes[0,0].set_title('(a) Scalar field')
axes[0,0].grid(True, alpha=0.3)

axes[0,1].plot(t, rho/RHO_C, 'r-')
axes[0,1].axhline(1.0, color='k', ls='--', label=r'$\rho_c$')
axes[0,1].set_xlabel('t'); axes[0,1].set_ylabel(r'$\rho/\rho_c$')
axes[0,1].set_title('(b) Density'); axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)

axes[1,0].plot(t, H, 'g-'); axes[1,0].axhline(0, color='k', ls=':')
axes[1,0].set_xlabel('t'); axes[1,0].set_ylabel('H')
axes[1,0].set_title('(c) Signed Hubble rate')
axes[1,0].grid(True, alpha=0.3)

axes[1,1].plot(t, a, 'm-')
axes[1,1].set_xlabel('t'); axes[1,1].set_ylabel('a')
axes[1,1].set_title('(d) Scale factor (bounce)')
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig10_lqc_bounce.png", dpi=100)
print("Saved: paper/figures/fig10_lqc_bounce.png")
