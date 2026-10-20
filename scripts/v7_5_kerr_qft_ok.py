"""v7.5-alpha2: Kerr + PTUM in the QFT-compatible window."""
import numpy as np, time
from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state
from zt005.bssn_refined_kerr_complex_ptum import evolve_refined_kerr_complex_ptum
from zt005.matter.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2); dt = 2e-4; steps = 100
print(f"Kerr + PTUM QFT-OK window")

for a in [0.0, 0.9]:
    for A in [1.6, 1.8]:
        t0 = time.time()
        p_c = kerr_complex_ptum_state(N=N_c, L=L_c, M=1.0, a=a,
                                        A_fuel=A, sigma_fuel=4.0)
        p_f = kerr_complex_ptum_state(N=N_f, L=L_f, M=1.0, a=a,
                                        A_fuel=A, sigma_fuel=4.0)
        params = p_c["ptum_params"]
        ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
        Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))

        sc = p_c["state"]; PRc = p_c["PR"]; PIc = p_c["PI"]
        Phc = p_c["Phi"]; psc = p_c["psi"]
        sf = p_f["state"]; PRf = p_f["PR"]; PIf = p_f["PI"]
        Phf = p_f["Phi"]; psf = p_f["psi"]

        t_hist = []; rho_hist = []
        for chunk_i in range(steps // 25):
            out = evolve_refined_kerr_complex_ptum(
                sc, PRc, PIc, Zc, Zc, Phc, Zc, psc, Zc,
                sf, PRf, PIf, Zf, Zf, Phf, Zf, psf, Zf,
                h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=25,
                ptum_params=params)
            sc = out["state_c"]; Phc = out["Ph_c"]; psc = out["ps_c"]
            sf = out["state_f"]; Phf = out["Ph_f"]; psf = out["ps_f"]
            step = (chunk_i+1)*25
            src = ptum_sources(Phc, Zc, psc, Zc, sc.gtilde, sc.phi,
                               h_c, **ptum_kwargs)
            t_hist.append(step*dt)
            rho_hist.append(float(src["rho"].min()))

        t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
        rho_min = float(np.min(rho_arr))
        r = check_ford_roman(rho_arr, t_arr, tau=1.0, N_field_components=2)
        print(f"a={a:.1f} A={A:.2f}: rho_min={rho_min:.4e} "
              f"ratio={r['ratio']:.4f} viol={r['violates']} "
              f"t={time.time()-t0:.0f}s")
