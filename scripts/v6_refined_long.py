import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_ptum import evolve_refined_ptum
from zt005.ptum_field import fuel_profile, psi_equilibrium
from zt005.matter_ptum import ptum_sources

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; total_steps = 1000; chunk = 50
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
A = 3.0

s_c = minkowski_state(N_c)
Phi_c = fuel_profile(N_c, L_c, Phi0=A, sigma=4.0)
psi_c = psi_equilibrium(Phi_c, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Zc = np.zeros((N_c, N_c, N_c))

s_f = minkowski_state(N_f)
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
Rf2 = Xf**2 + Yf**2 + Zf**2
Phi_f = A * np.exp(-Rf2 / (2.0 * 4.0**2))
psi_f = psi_equilibrium(Phi_f, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Zf = np.zeros((N_f, N_f, N_f))

print(f"REFINED LONG PTUM: {total_steps} steps, A={A}")
t0 = time.time()
n_neg_series = []
for chunk_i in range(total_steps // chunk):
    out = evolve_refined_ptum(
        s_c, Phi_c, Zc, psi_c, Zc,
        s_f, Phi_f, Zf, psi_f, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    s_c = out["state_c"]; Phi_c = out["Phi_c"]; psi_c = out["psi_c"]
    s_f = out["state_f"]; Phi_f = out["Phi_f"]; psi_f = out["psi_f"]
    PiPhi_c = out["PiPhi_c"]; Pipsi_c = out["Pipsi_c"]

    src = ptum_sources(Phi_c, PiPhi_c, psi_c, Pipsi_c, s_c.gtilde, s_c.phi,
                       h_c, V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                       lam_psi=lam_psi)
    rho_min = float(src["rho"].min())
    n_neg = int(np.sum(src["rho"] < 0))
    step = (chunk_i+1)*chunk
    n_neg_series.append((step, rho_min, n_neg))
    fin = (np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
           and np.all(np.isfinite(Phi_c)) and np.all(np.isfinite(Phi_f)))
    if step % 100 == 0:
        print(f"step={step:4d} rho_min={rho_min:.4e} neg={n_neg:3d} "
              f"Phimax={Phi_c.max():.3e} fin={fin} t={time.time()-t0:.0f}s")
    if not fin:
        print("NON-FINITE — stopping")
        break

# Save time series for Ford-Roman test
np.savez("results/validated/v6_long_ptum.npz",
         steps=np.array([s for s,_,_ in n_neg_series]),
         rho_min=np.array([r for _,r,_ in n_neg_series]),
         n_neg=np.array([n for _,_,n in n_neg_series]))
print("Saved: results/validated/v6_long_ptum.npz")
