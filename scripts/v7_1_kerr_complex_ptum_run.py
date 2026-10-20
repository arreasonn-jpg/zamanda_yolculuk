import numpy as np, time
from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state
from zt005.bssn_refined_kerr_complex_ptum import evolve_refined_kerr_complex_ptum
from zt005.matter.matter_ptum import ptum_sources

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 100; chunk = 25
A_fuel = 3.0; a = 0.9
print(f"REFINED KERR+COMPLEX+PTUM a={a} A={A_fuel}")

p_c = kerr_complex_ptum_state(N=N_c, L=L_c, M=1.0, a=a, A_fuel=A_fuel,
                                sigma_fuel=4.0)
p_f = kerr_complex_ptum_state(N=N_f, L=L_f, M=1.0, a=a, A_fuel=A_fuel,
                                sigma_fuel=4.0)
params = p_c["ptum_params"]
Zc = np.zeros((N_c,N_c,N_c))

ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
src0 = ptum_sources(p_c["Phi"], Zc, p_c["psi"], Zc,
                    p_c["state"].gtilde, p_c["state"].phi, h_c, **ptum_kwargs)
print(f"initial rho_min={src0['rho'].min():.3e}")

sc = p_c["state"]; PRc = p_c["PR"]; PIc = p_c["PI"]
Phc = p_c["Phi"]; psc = p_c["psi"]
sf = p_f["state"]; PRf = p_f["PR"]; PIf = p_f["PI"]
Phf = p_f["Phi"]; psf = p_f["psi"]
Zf = np.zeros((N_f,N_f,N_f))

t0 = time.time()
for chunk_i in range(steps // chunk):
    out = evolve_refined_kerr_complex_ptum(
        sc, PRc, PIc, Zc, Zc, Phc, Zc, psc, Zc,
        sf, PRf, PIf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        ptum_params=params)
    sc = out["state_c"]; PRc = out["PR_c"]; PIc = out["PI_c"]
    Phc = out["Ph_c"]; psc = out["ps_c"]
    sf = out["state_f"]; PRf = out["PR_f"]; PIf = out["PI_f"]
    Phf = out["Ph_f"]; psf = out["ps_f"]

    step = (chunk_i + 1) * chunk
    src = ptum_sources(Phc, Zc, psc, Zc, sc.gtilde, sc.phi, h_c, **ptum_kwargs)
    rho_min = float(src["rho"].min())
    n_neg = int(np.sum(src["rho"] < 0))
    Kmax_c = float(np.max(np.abs(sc.K)))
    Kmax_f = float(np.max(np.abs(sf.K)))
    fin = (np.all(np.isfinite(sc.phi)) and np.all(np.isfinite(sf.phi))
           and np.all(np.isfinite(Phc)) and np.all(np.isfinite(psc)))
    print(f"step={step:4d} rho_min={rho_min:.4e} n_neg={n_neg} "
          f"Kmax_c={Kmax_c:.4e} Kmax_f={Kmax_f:.4e} "
          f"Phmax={np.max(np.abs(Phc)):.4f} psimax={np.max(np.abs(psc)):.4f} "
          f"alpha_min={sc.alpha.min():.4f} fin={fin} t={time.time()-t0:.0f}s")
    if not fin: break
