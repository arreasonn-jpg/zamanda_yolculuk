"""C1: refined grid + lambda=1.0 fast discharge."""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_complex_ptum_v2 import evolve_refined_complex_ptum_v2
from zt005.matter.ptum_field import fuel_profile, psi_equilibrium_smooth
from zt005.matter.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
offset = (2, 2, 2)
dt = 2e-4; steps = 100; chunk = 50
A_fuel = 1.5
params = {"V0": 1.0, "lam": 1.0, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

print(f"C1 REFINED FAST DISCHARGE: A={A_fuel} lambda=1.0")

xc = np.linspace(-L_c/2, L_c/2, N_c)
xf = np.linspace(-L_f/2, L_f/2, N_f)
Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing='ij')
Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing='ij')
Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2

s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
Phc = A_fuel*np.exp(-Rc2/8.0); Phf = A_fuel*np.exp(-Rf2/8.0)
psc = psi_equilibrium_smooth(Phc); psf = psi_equilibrium_smooth(Phf)
Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))

t_hist=[]; rho_hist=[]; phimax_hist=[]
t0 = time.time()
for chunk_i in range(steps // chunk):
    out = evolve_refined_complex_ptum_v2(
        s_c, Zc, Zc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, Zf, Zf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=offset, dt=dt, n_steps=chunk,
        ptum_params=params, restrict_matter=False)
    s_c = out["state_c"]; Phc = out["Ph_c"]; psc = out["ps_c"]
    s_f = out["state_f"]; Phf = out["Ph_f"]; psf = out["ps_f"]
    step = (chunk_i+1)*chunk
    src = ptum_sources(Phc, Zc, psc, Zc, s_c.gtilde, s_c.phi, h_c, **ptum_kwargs)
    t_hist.append(step*dt)
    rho_hist.append(float(src["rho"].min()))
    phimax_hist.append(float(np.max(np.abs(Phc))))
    if step % 20 == 0:
        print(f"step={step:4d} t={step*dt:.3f} rho_min={rho_hist[-1]:.4e} "
              f"Phimax={phimax_hist[-1]:.4f} psimax={np.max(np.abs(psc)):.4f} "
              f"Kmax_c={np.max(np.abs(s_c.K)):.3e} t_wall={time.time()-t0:.0f}s")

np.savez("results/validated/v7_2_c1_refined_fast.npz",
         t=np.array(t_hist), rho=np.array(rho_hist),
         phimax=np.array(phimax_hist))

t_arr=np.array(t_hist); rho_arr=np.array(rho_hist)
print("\n=== Ford-Roman massless (N=2) ===")
for tau in [0.1, 0.5, 1.0]:
    r = check_ford_roman(rho_arr, t_arr, tau=tau, N_field_components=2)
    print(f"tau={tau:.2f} value={r['value']:.4e} bound={r['bound']:.4e} "
          f"ratio={r['ratio']:.4f} violates={r['violates']}")

# Phi below 1
below = np.array(phimax_hist) < 1.0
if below.any():
    idx = int(np.argmax(below))
    print(f"\nPhi < 1 at t = {t_hist[idx]:.3f}")
