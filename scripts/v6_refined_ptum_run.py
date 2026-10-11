import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_ptum import evolve_refined_ptum
from zt005.ptum_field import fuel_profile, m_eff_squared, psi_equilibrium
from zt005.matter_ptum import ptum_sources

# Coarse: wide domain, coarse grid
N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
# Fine: centered, 2:1 refinement
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)  # coarse index of fine cell (0,0,0)

dt = 2e-4; steps = 200
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
A_fuel = 3.0  # the amplitude that failed in v5.0

print(f"REFINED PTUM RUN: N_c={N_c} h_c={h_c:.2f} N_f={N_f} h_f={h_f:.2f}")
print(f"dt={dt} steps={steps} A_fuel={A_fuel}")

# Coarse initial data
s_c = minkowski_state(N_c)
Phi_c = fuel_profile(N_c, L_c, Phi0=A_fuel, sigma=4.0)
psi_c = psi_equilibrium(Phi_c, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Zc = np.zeros((N_c, N_c, N_c))

# Fine initial data (same profile, finer sampling)
s_f = minkowski_state(N_f)
xc = np.linspace(-L_c/2, L_c/2, N_c)
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
Rf = Xf**2 + Yf**2 + Zf**2
Phi_f = A_fuel * np.exp(-Rf / (2.0 * 4.0**2))
psi_f = psi_equilibrium(Phi_f, m2=m2, alpha=alpha_p, Eth=Eth, lam_psi=lam_psi)
Zf = np.zeros((N_f, N_f, N_f))

# Initial diagnostics
meff_c = m_eff_squared(Phi_c, m2, alpha_p, Eth)
meff_f = m_eff_squared(Phi_f, m2, alpha_p, Eth)
print(f"Activated coarse: {np.sum(meff_c < 0)}/{Phi_c.size}")
print(f"Activated fine:   {np.sum(meff_f < 0)}/{Phi_f.size}")

src_c = ptum_sources(Phi_c, Zc, psi_c, Zc, s_c.gtilde, s_c.phi, h_c,
                     V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                     lam_psi=lam_psi)
print(f"rho_min coarse: {src_c['rho'].min():.3e}")

# Evolution: run in chunks to inspect progress
chunk = 25
t0 = time.time()
for chunk_i in range(steps // chunk):
    out = evolve_refined_ptum(
        s_c, Phi_c, Zc, psi_c, Zc,
        s_f, Phi_f, Zf, psi_f, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    s_c = out["state_c"]; Phi_c = out["Phi_c"]; psi_c = out["psi_c"]
    s_f = out["state_f"]; Phi_f = out["Phi_f"]; psi_f = out["psi_f"]

    src_c = ptum_sources(Phi_c, Zc, psi_c, Zc, s_c.gtilde, s_c.phi, h_c,
                         V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                         lam_psi=lam_psi)
    rho_min_c = src_c['rho'].min()
    Kmax_c = np.max(np.abs(s_c.K))
    Kmax_f = np.max(np.abs(s_f.K))
    fin = (np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
           and np.all(np.isfinite(Phi_c)) and np.all(np.isfinite(Phi_f)))
    print(f"step={(chunk_i+1)*chunk:4d} rho_min_c={rho_min_c:.3e} "
          f"Kmax_c={Kmax_c:.3e} Kmax_f={Kmax_f:.3e} finite={fin} "
          f"t={time.time()-t0:.0f}s")
    if not fin:
        print("FAILED — non-finite detected")
        break
