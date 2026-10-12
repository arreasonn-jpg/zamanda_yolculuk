import numpy as np, time
from zt005.kerr_sc_ptum_ic import kerr_sc_ptum_state, combined_diagnostics
from zt005.bssn_rk4_sc_ptum import rk4_step_sc_ptum

N = 9; L = 20.0; h = L/(N-1); dt = 1e-4; steps = 30
print(f"FULL SC+PTUM EVOLUTION: N={N} dt={dt} steps={steps}")

for a in [0.0, 0.5, 0.9]:
    pack = kerr_sc_ptum_state(N=N, L=L, M=1.0, a=a, rc=2.0,
                                A_fuel=1.5, sigma_fuel=4.0)
    s = pack["state"]
    phi_sc = pack["phi_sc"]; Pi_sc = pack["Pi_sc"]
    Phi = pack["Phi"]; Pi_Phi = pack["Pi_Phi"]
    psi = pack["psi"]; Pi_psi = pack["Pi_psi"]
    V_sc = pack["V_sc"]
    ptum_params = pack["V_ptum_params"]

    t0 = time.time()
    fin = True
    for step in range(steps):
        try:
            s, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi = rk4_step_sc_ptum(
                s, phi_sc, Pi_sc, Phi, Pi_Phi, psi, Pi_psi,
                h, dt, V_sc, ptum_params)
        except Exception as e:
            fin = False
            print(f"  a={a} step {step} ERROR: {e}")
            break

    pack["state"] = s
    pack["phi_sc"] = phi_sc; pack["Pi_sc"] = Pi_sc
    pack["Phi"] = Phi; pack["Pi_Phi"] = Pi_Phi
    pack["psi"] = psi; pack["Pi_psi"] = Pi_psi

    d = combined_diagnostics(pack, h)
    fin = fin and (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
                    and np.all(np.isfinite(phi_sc)))
    print(f"a={a:.2f} rho_min={d['rho_total_min']:.3e} "
          f"rho_sc={d['rho_sc_min']:.3e} rho_pt={d['rho_ptum_min']:.3e} "
          f"Kmax={d['Kmax']:.3e} alpha_min={d['alpha_min']:.4f} "
          f"fin={fin} t={time.time()-t0:.0f}s")
