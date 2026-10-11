import numpy as np, time
from zt005.kerr_sc_ic import kerr_sc_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.bssn_ricci_full import hamiltonian_constraint_full

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 50
print(f"KERR + SC CORE: N={N} dt={dt} steps={steps}")
print(f"{'a':>5} {'rc':>5} {'H0_ext':>11} {'Hf_ext':>11} {'growth':>8} "
      f"{'alpha_min':>10} {'finite':>7}")

# mask for exterior (outside core)
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)

for rc in [1.0, 2.0]:
    for a in [0.0, 0.5, 0.9]:
        s = kerr_sc_state(N=N, L=L, M=1.0, a=a, rc=rc)
        # Zero matter field
        phi_f = np.zeros((N,N,N))
        Pi_f = np.zeros((N,N,N))
        mask = R > 3.0*rc
        H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        t0 = time.time()
        for step in range(steps):
            s, phi_f, Pi_f = rk4_step_matter(s, phi_f, Pi_f, h, dt)
        Hf = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.gtilde)))
        print(f"{a:5.2f} {rc:5.2f} {H0:>11.4e} {Hf:>11.4e} {Hf/H0:>8.4f} "
              f"{s.alpha.min():>10.4f} {str(fin):>7}")
