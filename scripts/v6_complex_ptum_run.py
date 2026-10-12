import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum

N = 9; L = 20.0; h = L/(N-1); dt = 1e-4; steps = 30
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
print(f"COMPLEX + PTUM: N={N} dt={dt} steps={steps}")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

for A in [0.5, 1.5]:
    s = minkowski_state(N)
    PR = 1e-4 * np.exp(-R2/2.0)
    PI = 0.5e-4 * np.exp(-R2/2.0)
    PiR = np.zeros((N,N,N))
    PiI = 0.5e-4 * np.exp(-R2/2.0)
    Ph = A * np.exp(-R2/8.0)
    PPh = np.zeros((N,N,N))
    ps = np.zeros((N,N,N))
    Pps = np.zeros((N,N,N))

    t0 = time.time()
    fin = True
    for step in range(steps):
        s, PR, PI, PiR, PiI, Ph, PPh, ps, Pps = rk4_step_complex_ptum(
            s, PR, PI, PiR, PiI, Ph, PPh, ps, Pps,
            h, dt, params)
    Kmax = float(np.max(np.abs(s.K)))
    phimax = float(np.max(np.abs(Ph)))
    psimax = float(np.max(np.abs(ps)))
    fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(PR))
           and np.all(np.isfinite(Ph)) and np.all(np.isfinite(ps)))
    print(f"A={A:.2f} Kmax={Kmax:.4e} Phimax={phimax:.4e} "
          f"psimax={psimax:.4e} fin={fin} t={time.time()-t0:.0f}s")
