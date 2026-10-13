"""v8.0-alpha4: does Brans-Dicke phi shift the PTUM QFT-OK window?"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_bd_ptum import rk4_step_bd_ptum
from zt005.ptum_field import psi_equilibrium_smooth
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
omega = 40000.0  # Cassini bound
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

print(f"BD + PTUM scan: omega={omega}")
print(f"{'phi_BD':>8} {'A':>5} {'rho_min':>12} {'ratio':>10} {'viol':>6}")

# phi_BD modifies G_eff = G/phi
for phi_BD_val in [0.5, 1.0, 2.0]:
    for A in [1.6, 1.8, 2.0]:
        t0 = time.time()
        s = minkowski_state(N)
        Z = np.zeros((N,N,N))
        phi_BD = np.ones((N,N,N)) * phi_BD_val
        Pi_BD = np.zeros((N,N,N))
        Phi = A*np.exp(-R2/8.0)
        psi = psi_equilibrium_smooth(Phi)

        t_hist = []; rho_hist = []
        for step in range(steps):
            s, phi_BD, Pi_BD, Phi, Z, psi, Z = rk4_step_bd_ptum(
                s, phi_BD, Pi_BD, Phi, Z, psi, Z, h, dt, omega, params)
            if (step + 1) % 5 == 0:
                src = ptum_sources(Phi, Z, psi, Z, s.gtilde, s.phi, h, **ptum_kwargs)
                t_hist.append((step+1)*dt)
                rho_hist.append(float(src["rho"].min()))
        t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
        rho_min = float(np.min(rho_arr))
        r = check_ford_roman(rho_arr, t_arr, tau=1.0, N_field_components=2)
        print(f"{phi_BD_val:>8.2f} {A:>5.2f} {rho_min:>12.4e} {r['ratio']:>10.4f} "
              f"{str(r['violates']):>6} t={time.time()-t0:.0f}s")
