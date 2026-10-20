import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_sc_ptum import evolve_refined_sc_ptum
from zt005.matter.sc_condensate import sc_condensate_potential
from zt005.matter.matter_scalar import scalar_sources
from zt005.matter.matter_ptum import ptum_sources
from zt005.kerr_sc_ptum_ic import kerr_sc_ptum_state

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 100; chunk = 25
a = 0.9
A_fuel = 3.0  # the strong amplitude

print(f"REFINED FULL STACK (STRONG): A_fuel={A_fuel} N_c={N_c} N_f={N_f}")

pk_c = kerr_sc_ptum_state(N=N_c, L=L_c, M=1.0, a=a, rc=2.0,
                           A_fuel=A_fuel, sigma_fuel=4.0)
pk_f = kerr_sc_ptum_state(N=N_f, L=L_f, M=1.0, a=a, rc=2.0,
                           A_fuel=A_fuel, sigma_fuel=4.0)

s_c = pk_c["state"]; psc_c = pk_c["phi_sc"]; Psc_c = pk_c["Pi_sc"]
Ph_c = pk_c["Phi"]; PPh_c = pk_c["Pi_Phi"]; ps_c = pk_c["psi"]; Pps_c = pk_c["Pi_psi"]
s_f = pk_f["state"]; psc_f = pk_f["phi_sc"]; Psc_f = pk_f["Pi_sc"]
Ph_f = pk_f["Phi"]; PPh_f = pk_f["Pi_Phi"]; ps_f = pk_f["psi"]; Pps_f = pk_f["Pi_psi"]
V_sc = pk_c["V_sc"]
params = pk_c["V_ptum_params"]

# Diagnostics helper
def diag(s_c, psc_c, Psc_c, Ph_c, PPh_c, ps_c, Pps_c, h_c):
    src_sc = scalar_sources(psc_c, Psc_c, s_c.gtilde, s_c.phi, h_c, V_func=V_sc)
    ptum_kwargs = dict(params)
    if "alpha_p" in ptum_kwargs:
        ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
    src_pt = ptum_sources(Ph_c, PPh_c, ps_c, Pps_c, s_c.gtilde, s_c.phi, h_c,
                          **ptum_kwargs)
    rho_tot = src_sc["rho"] + src_pt["rho"]
    return rho_tot, src_pt["rho"]

rho_tot, rho_pt = diag(s_c, psc_c, Psc_c, Ph_c, PPh_c, ps_c, Pps_c, h_c)
print(f"initial rho_total_min={rho_tot.min():.3e} rho_ptum_min={rho_pt.min():.3e}")

t0 = time.time()
for chunk_i in range(steps // chunk):
    try:
        out = evolve_refined_sc_ptum(
            s_c, psc_c, Psc_c, Ph_c, PPh_c, ps_c, Pps_c,
            s_f, psc_f, Psc_f, Ph_f, PPh_f, ps_f, Pps_f,
            h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
            V_sc=V_sc, ptum_params=params)
    except Exception as e:
        print(f"ERROR at step {(chunk_i+1)*chunk}: {e}")
        break

    s_c = out["state_c"]; psc_c = out["phi_sc_c"]; Psc_c = out["Pi_sc_c"]
    Ph_c = out["Phi_c"]; PPh_c = out["Pi_Phi_c"]
    ps_c = out["psi_c"]; Pps_c = out["Pi_psi_c"]
    s_f = out["state_f"]

    step = (chunk_i + 1) * chunk
    rho_tot, rho_pt = diag(s_c, psc_c, Psc_c, Ph_c, PPh_c, ps_c, Pps_c, h_c)
    Kmax_c = float(np.max(np.abs(s_c.K)))
    Kmax_f = float(np.max(np.abs(s_f.K)))
    fin = (np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
           and np.all(np.isfinite(Ph_c)) and np.all(np.isfinite(psc_c)))
    n_neg = int(np.sum(rho_tot < 0))
    print(f"step={step:4d} rho_min={rho_tot.min():.4e} rho_pt={rho_pt.min():.4e} "
          f"n_neg={n_neg} Kmax_c={Kmax_c:.3e} Kmax_f={Kmax_f:.3e} "
          f"alpha_min={s_c.alpha.min():.4f} fin={fin} t={time.time()-t0:.0f}s")
    if not fin:
        print("NON-FINITE")
        break
