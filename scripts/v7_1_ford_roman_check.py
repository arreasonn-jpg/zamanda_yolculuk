"""Ford-Roman check on the v7.1 Kerr+complex+PTUM config."""
import numpy as np, time
from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state
from zt005.bssn_refined_kerr_complex_ptum import evolve_refined_kerr_complex_ptum
from zt005.matter.matter_ptum import ptum_sources
from zt005.quantum_inequality import (
    check_ford_roman, check_ford_roman_massive,
)

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 100; chunk = 25
A_fuel = 3.0; a = 0.9
print(f"Ford-Roman check: Kerr+complex+PTUM a={a} A={A_fuel}")

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

t_hist = []; rho_hist = []
t0 = time.time()
for chunk_i in range(steps // chunk):
    out = evolve_refined_kerr_complex_ptum(
        sc, PRc, PIc, Zc, Zc, Phc, Zc, psc, Zc,
        sf, PRf, PIf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        ptum_params=params)
    sc = out["state_c"]; Phc = out["Ph_c"]; psc = out["ps_c"]
    sf = out["state_f"]; Phf = out["Ph_f"]; psf = out["ps_f"]
    step = (chunk_i + 1) * chunk
    src = ptum_sources(Phc, Zc, psc, Zc, sc.gtilde, sc.phi, h_c, **ptum_kwargs)
    rho_hist.append(float(src["rho"].min()))
    t_hist.append(step * dt)
print(f"run done in {time.time()-t0:.0f}s")
print(f"rho_min: first={rho_hist[0]:.4e} last={rho_hist[-1]:.4e}")

t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
N_comp = 2  # PTUM has two real scalars

print("\n=== Massless Ford-Roman (N=2) ===")
print(f"{'tau':>8} {'value':>14} {'bound':>14} {'ratio':>10} {'violates':>10}")
for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
    r = check_ford_roman(rho_arr, t_arr, tau=tau, N_field_components=N_comp)
    print(f"{tau:>8.3f} {r['value']:>14.4e} {r['bound']:>14.4e} "
          f"{r['ratio']:>10.4f} {str(r['violates']):>10}")

print("\n=== Massive Ford-Roman (mu=1.0, N=2) ===")
print(f"{'tau':>8} {'mu*tau':>8} {'value':>14} {'bound':>14} {'ratio':>10} {'violates':>10}")
for tau in [0.01, 0.05, 0.1, 0.5, 1.0]:
    r = check_ford_roman_massive(rho_arr, t_arr, tau=tau, mu=1.0,
                                  N_field_components=N_comp)
    print(f"{tau:>8.3f} {r['mu_tau']:>8.3f} {r['value']:>14.4e} "
          f"{r['bound']:>14.4e} {r['ratio']:>10.4f} {str(r['violates']):>10}")
