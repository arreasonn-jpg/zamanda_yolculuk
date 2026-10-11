import numpy as np, time
from zt005.kerr_ptum_ic import kerr_ptum_state, kerr_ptum_diagnostics
from zt005.bssn_rk4_ptum import rk4_step_ptum

N = 9; L = 20.0; M = 1.0; h = L/(N-1); dt = 2e-4; steps = 50
A_fuel = 3.0; sigma = 4.0
V0=1.0; lam=0.05; m2=1.0; alpha_p=1.0; Eth=1.0; lam_psi=1.0
print(f"KERR + PTUM: N={N} dt={dt} steps={steps} A_fuel={A_fuel}")
print(f"{'a':>5} {'rho_min':>11} {'n_neg':>6} {'n_act':>6} "
      f"{'Kmax':>11} {'alpha_min':>10} {'finite':>7} {'t':>4}")

for a in [0.0, 0.5, 0.9]:
    s, Phi, psi, PiPhi, Pipsi = kerr_ptum_state(
        N=N, L=L, M=M, a=a, A_fuel=A_fuel, sigma=sigma,
        m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    t0 = time.time()
    for step in range(steps):
        s, Phi, PiPhi, psi, Pipsi = rk4_step_ptum(
            s, Phi, PiPhi, psi, Pipsi, h, dt,
            V0=V0, lam=lam, m2=m2, alpha_p=alpha_p, Eth=Eth, lam_psi=lam_psi)
    d = kerr_ptum_diagnostics(s, Phi, PiPhi, psi, Pipsi, h,
                               V0=V0, lam=lam, m2=m2, alpha_p=alpha_p,
                               Eth=Eth, lam_psi=lam_psi)
    fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(Phi))
           and np.all(np.isfinite(psi)))
    print(f"{a:5.2f} {d['rho_min']:>11.3e} {d['n_neg']:>6d} "
          f"{d['n_active']:>6d} {d['Kmax']:>11.3e} "
          f"{d['alpha_min']:>10.4f} {str(fin):>7} {time.time()-t0:>4.0f}")
