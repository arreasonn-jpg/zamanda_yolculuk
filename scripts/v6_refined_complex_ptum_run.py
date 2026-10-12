import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_complex_ptum import evolve_refined_complex_ptum
from zt005.ptum_field import psi_equilibrium
from zt005.matter_ptum import ptum_sources

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 100; chunk = 25
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
A_fuel = 3.0
print(f"REFINED COMPLEX+PTUM A={A_fuel}")

xc = np.linspace(-L_c/2, L_c/2, N_c)
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing="ij")
Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing="ij")
Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2

s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
PRc = 1e-4*np.exp(-Rc2/2.0); PRf = 1e-4*np.exp(-Rf2/2.0)
PIc = 0.5e-4*np.exp(-Rc2/2.0); PIf = 0.5e-4*np.exp(-Rf2/2.0)
Phc = A_fuel*np.exp(-Rc2/8.0); Phf = A_fuel*np.exp(-Rf2/2.0)
psc = psi_equilibrium(Phc); psf = psi_equilibrium(Phf)
Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))

t0 = time.time()
for chunk_i in range(steps // chunk):
    out = evolve_refined_complex_ptum(
        s_c, PRc, PIc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, PRf, PIf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        ptum_params=params)
    s_c = out["state_c"]; PRc = out["PR_c"]; PIc = out["PI_c"]
    Phc = out["Ph_c"]; psc = out["ps_c"]
    s_f = out["state_f"]; PRf = out["PR_f"]; PIf = out["PI_f"]
    Phf = out["Ph_f"]; psf = out["ps_f"]

    step = (chunk_i + 1) * chunk
    ptum_kwargs = dict(params)
    ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")
    src = ptum_sources(Phc, Zc, psc, Zc, s_c.gtilde, s_c.phi, h_c, **ptum_kwargs)
    rho_min = float(src["rho"].min())
    n_neg = int(np.sum(src["rho"] < 0))
    Kmax_c = float(np.max(np.abs(s_c.K)))
    Kmax_f = float(np.max(np.abs(s_f.K)))
    fin = (np.all(np.isfinite(s_c.phi)) and np.all(np.isfinite(s_f.phi))
           and np.all(np.isfinite(Phc)) and np.all(np.isfinite(psc)))
    print(f"step={step:4d} rho_min={rho_min:.4e} n_neg={n_neg} "
          f"Kmax_c={Kmax_c:.4e} Kmax_f={Kmax_f:.4e} "
          f"alpha_min={s_c.alpha.min():.4f} fin={fin} t={time.time()-t0:.0f}s")
    if not fin:
        print("NON-FINITE")
        break
