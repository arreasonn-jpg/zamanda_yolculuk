import numpy as np, time
from zt005.refined_kerr_ptum_ic import refined_kerr_ptum_ic
from zt005.bssn_refined_ptum import evolve_refined_ptum
from zt005.matter.matter_ptum import ptum_sources

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; total_steps = 500; chunk = 50
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
A_fuel = 3.0; sigma = 4.0

print(f"REFINED KERR+PTUM LONG: {total_steps} steps, a=0.9")
a = 0.9

s_c, Phi_c, psi_c, s_f, Phi_f, psi_f = refined_kerr_ptum_ic(
    N_c=N_c, L_c=L_c, N_f=N_f, L_f=L_f, M=1.0, a=a, A_fuel=A_fuel,
    sigma_phys=sigma, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
Zc = np.zeros((N_c, N_c, N_c)); Zf = np.zeros((N_f, N_f, N_f))

t0 = time.time()
for chunk_i in range(total_steps // chunk):
    out = evolve_refined_ptum(
        s_c, Phi_c, Zc, psi_c, Zc, s_f, Phi_f, Zf, psi_f, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    s_c = out["state_c"]; Phi_c = out["Phi_c"]; psi_c = out["psi_c"]
    s_f = out["state_f"]; Phi_f = out["Phi_f"]; psi_f = out["psi_f"]
    PiPhi_c = out["PiPhi_c"]; Pipsi_c = out["Pipsi_c"]

    step = (chunk_i + 1) * chunk
    src = ptum_sources(Phi_c, PiPhi_c, psi_c, Pipsi_c, s_c.gtilde, s_c.phi,
                       h_c, V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                       lam_psi=lam_psi)
    rho_min = float(src["rho"].min())
    n_neg = int(np.sum(src["rho"] < 0))
    Kmax_c = float(np.max(np.abs(s_c.K)))
    Kmax_f = float(np.max(np.abs(s_f.K)))
    fin = (np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
           and np.all(np.isfinite(Phi_c)) and np.all(np.isfinite(Phi_f)))
    if step % 100 == 0:
        print(f"step={step:4d} rho_min={rho_min:.4e} neg={n_neg:3d} "
              f"Kmax_c={Kmax_c:.3e} Kmax_f={Kmax_f:.3e} fin={fin} "
              f"t={time.time()-t0:.0f}s")
    if not fin:
        print("NON-FINITE — stopping")
        break
