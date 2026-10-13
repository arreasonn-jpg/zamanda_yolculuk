"""v7.2: transient PTUM discharge with natural Phi decay.

Uses psi_from_equilibrium=False so Phi decays by its own
potential. Measures the discharge timescale and applies Ford-Roman
over the full time series (which now includes both discharge and
post-discharge).
"""
import numpy as np, time
from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state
from zt005.bssn_refined_kerr_complex_ptum import evolve_refined_kerr_complex_ptum
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman, check_ford_roman_massive

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 500; chunk = 50
A_fuel = 3.0; a = 0.9
print(f"v7.2 TRANSIENT DISCHARGE: A={A_fuel} a={a} steps={steps}")

p_c = kerr_complex_ptum_state(N=N_c, L=L_c, M=1.0, a=a, A_fuel=A_fuel,
                                sigma_fuel=4.0)
p_f = kerr_complex_ptum_state(N=N_f, L=L_f, M=1.0, a=a, A_fuel=A_fuel,
                                sigma_fuel=4.0)
params = p_c["ptum_params"]
Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

sc = p_c["state"]; PRc = p_c["PR"]; PIc = p_c["PI"]
Phc = p_c["Phi"]; psc = p_c["psi"]
sf = p_f["state"]; PRf = p_f["PR"]; PIf = p_f["PI"]
Phf = p_f["Phi"]; psf = p_f["psi"]

t_hist = []; rho_hist = []; phimax_hist = []
t0 = time.time()
for chunk_i in range(steps // chunk):
    # NOTE: psi_from_equilibrium=False -> natural Phi decay
    out = evolve_refined_kerr_complex_ptum(
        sc, PRc, PIc, Zc, Zc, Phc, Zc, psc, Zc,
        sf, PRf, PIf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        ptum_params=params, psi_from_equilibrium=False)
    sc = out["state_c"]; Phc = out["Ph_c"]; psc = out["ps_c"]
    sf = out["state_f"]; Phf = out["Ph_f"]; psf = out["ps_f"]
    step = (chunk_i + 1) * chunk
    src = ptum_sources(Phc, Zc, psc, Zc, sc.gtilde, sc.phi, h_c, **ptum_kwargs)
    rho_hist.append(float(src["rho"].min()))
    t_hist.append(step * dt)
    phimax_hist.append(float(np.max(np.abs(Phc))))
    print(f"step={step:4d} t={step*dt:.4f} rho_min={rho_hist[-1]:.4e} "
          f"Phimax={phimax_hist[-1]:.4f} psimax={np.max(np.abs(psc)):.4f} "
          f"Kmax_c={np.max(np.abs(sc.K)):.3e} t_wall={time.time()-t0:.0f}s")

# Save time series
np.savez("results/validated/v7_2_transient_discharge.npz",
         t=np.array(t_hist), rho=np.array(rho_hist),
         phimax=np.array(phimax_hist))

# Ford-Roman test
t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
N_comp = 2
print("\n=== Massless Ford-Roman ===")
print(f"{'tau':>8} {'value':>14} {'bound':>14} {'ratio':>10} {'viol':>6}")
for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
    r = check_ford_roman(rho_arr, t_arr, tau=tau, N_field_components=N_comp)
    print(f"{tau:>8.3f} {r['value']:>14.4e} {r['bound']:>14.4e} "
          f"{r['ratio']:>10.4f} {str(r['violates']):>6}")

print("\n=== Massive (mu=1) Ford-Roman ===")
print(f"{'tau':>8} {'mu*tau':>8} {'value':>14} {'bound':>14} {'ratio':>10} {'viol':>6}")
for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
    r = check_ford_roman_massive(rho_arr, t_arr, tau=tau, mu=1.0,
                                  N_field_components=N_comp)
    print(f"{tau:>8.3f} {r['mu_tau']:>8.3f} {r['value']:>14.4e} "
          f"{r['bound']:>14.4e} {r['ratio']:>10.4f} {str(r['violates']):>6}")

# Discharge duration: how long until rho > -1?
threshold = -1.0
below = rho_arr < threshold
first_below = np.argmax(below)
last_below = len(below) - 1 - np.argmax(below[::-1])
print(f"\nDischarge (rho < {threshold}) from t={t_arr[first_below]:.4f} "
      f"to t={t_arr[last_below]:.4f}, duration={t_arr[last_below]-t_arr[first_below]:.4f}")
