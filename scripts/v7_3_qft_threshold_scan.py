"""v7.3: amplitude scan to find QFT-compatibility threshold.

For each A_fuel, run 100 steps of single-grid complex+PTUM,
extract rho(t), apply Ford-Roman (massless, N=2), report
tau=1.0 ratio.

Goal: find rho_min threshold where QFT compatibility breaks.
"""
import numpy as np, time
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex_ptum import rk4_step_complex_ptum
from zt005.ptum_field import fuel_profile, psi_equilibrium_smooth
from zt005.matter_ptum import ptum_sources
from zt005.quantum_inequality import check_ford_roman

N = 9; L = 20.0; h = L/(N-1); dt = 5e-4; steps = 100
params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
          "Eth": 1.0, "lam_psi": 1.0}
ptum_kwargs = dict(params); ptum_kwargs["alpha"] = ptum_kwargs.pop("alpha_p")

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z2**2

print(f"QFT threshold scan: N={N} dt={dt} steps={steps}")
print(f"{'A':>5} {'rho_min':>12} {'ratio_tau1':>12} {'violates':>10} {'t_wall':>8}")

results = []
for A in [1.5, 2.0, 2.5, 3.0, 3.5]:
    t0 = time.time()
    s = minkowski_state(N)
    Z = np.zeros((N,N,N))
    Phi = A*np.exp(-R2/8.0)
    psi = psi_equilibrium_smooth(Phi)
    PR = np.zeros((N,N,N)); PI = np.zeros((N,N,N))
    # Small complex scalar seed (mimics v7.1)
    PR = 1e-4*np.exp(-R2/2.0)
    PI = 0.5e-4*np.exp(-R2/2.0)

    t_hist = []; rho_hist = []
    try:
        for step in range(steps):
            s, PR, PI, Z, Z, Phi, Z, psi, Z = rk4_step_complex_ptum(
                s, PR, PI, Z, Z, Phi, Z, psi, Z, h, dt, params)
            if (step + 1) % 5 == 0:
                src = ptum_sources(Phi, Z, psi, Z, s.gtilde, s.phi, h, **ptum_kwargs)
                t_hist.append((step+1)*dt)
                rho_hist.append(float(src["rho"].min()))
        t_arr = np.array(t_hist); rho_arr = np.array(rho_hist)
        rho_min = float(np.min(rho_arr))
        r = check_ford_roman(rho_arr, t_arr, tau=1.0, N_field_components=2)
        ratio = r["ratio"]
        viol = r["violates"]
        print(f"{A:>5.2f} {rho_min:>12.4e} {ratio:>12.4f} {str(viol):>10} {time.time()-t0:>8.0f}")
        results.append((A, rho_min, ratio, viol))
    except Exception as e:
        print(f"{A:>5.2f} ERROR: {e}")
        break

# Analysis
print("\n=== Threshold analysis ===")
for A, rho_min, ratio, viol in results:
    marker = " <-- violation" if viol else ""
    print(f"A={A:.2f}: rho_min={rho_min:.3e}, ratio={ratio:.4f}{marker}")

# Find first violating A
first_viol = next((A for A, _, _, v in results if v), None)
last_ok = next((A for A, _, _, v in reversed(results) if not v), None)
if first_viol and last_ok:
    print(f"\nQFT compatibility threshold: {last_ok} < A_threshold < {first_viol}")
