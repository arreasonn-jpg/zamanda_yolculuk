"""v9.3-alpha1: quadrupole GW source.

Two scalar Gaussian bumps offset along x-axis (binary-like).
Non-axisymmetric source -> time-varying quadrupole -> GW emission.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.gw_extraction import (
    extract_at_observer, psi4_from_history,
    strain_amplitude_from_psi4_history,
)

N = 11; L = 20.0; h = L/(N-1); dt = 5e-3; steps = 100
A_pert = 1e-3; sigma = 1.5; d_sep = 3.0

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

# Two Gaussians offset along x
R2a = (X - d_sep/2)**2 + Y**2 + Z**2
R2b = (X + d_sep/2)**2 + Y**2 + Z**2
phi_sc = A_pert * (np.exp(-R2a/(2*sigma**2)) + np.exp(-R2b/(2*sigma**2)))
Pi_sc = np.zeros((N,N,N))
V_func = lambda p: 0.5 * p**2

s = minkowski_state(N)
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])

# Observer along z (perpendicular to binary axis)
obs_list = [(5, 5, 8), (5, 5, 9), (5, 5, 10)]

hp_hist = {obs: [] for obs in obs_list}
hc_hist = {obs: [] for obs in obs_list}
t_hist = []

print(f"Quadrupole GW: N={N}, steps={steps}, sep={d_sep}")

for step in range(steps):
    s, phi_sc, Pi_sc = rk4_step_matter(s, phi_sc, Pi_sc, h, dt, V_func=V_func)
    if step % 2 == 0:
        t_hist.append((step+1)*dt)
        for obs in obs_list:
            hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
            hp_hist[obs].append(hp)
            hc_hist[obs].append(hc)

t_arr = np.array(t_hist)
dt_out = t_arr[1] - t_arr[0]

results = {}
for obs in obs_list:
    hp = np.array(hp_hist[obs]); hc = np.array(hc_hist[obs])
    psi4 = psi4_from_history(hp, hc, dt_out)
    h_est = strain_amplitude_from_psi4_history(psi4, dt_out)
    results[obs] = {"hp": hp, "hc": hc, "psi4": psi4, "h_est": h_est}
    print(f"obs={obs}: h+ max={np.max(np.abs(hp)):.3e} "
          f"hx max={np.max(np.abs(hc)):.3e} "
          f"psi4 max={np.max(np.abs(psi4)) if len(psi4)>0 else 0:.3e}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
obs0 = obs_list[-1]
r = results[obs0]

axes[0,0].plot(t_arr, r["hp"], 'b-', label='h+')
axes[0,0].plot(t_arr, r["hc"], 'r--', label='hx')
axes[0,0].set_xlabel('t'); axes[0,0].set_ylabel('h')
axes[0,0].set_title(f'Strain at observer {obs0}')
axes[0,0].legend(); axes[0,0].grid(True, alpha=0.3)

if len(r["psi4"]) > 0:
    t_psi = t_arr[:len(r["psi4"])]
    axes[0,1].plot(t_psi, np.abs(r["psi4"]), 'g-')
    axes[0,1].set_xlabel('t'); axes[0,1].set_ylabel('|psi_4|')
    axes[0,1].set_title('NP scalar psi_4')
    axes[0,1].grid(True, alpha=0.3)

if len(r["h_est"]) > 0:
    t_h = t_arr[:len(r["h_est"])]
    axes[1,0].plot(t_h, r["h_est"], 'm-')
    axes[1,0].set_xlabel('t'); axes[1,0].set_ylabel('|h| recovered')
    axes[1,0].set_title('Strain recovery from psi_4')
    axes[1,0].grid(True, alpha=0.3)

for obs in obs_list:
    axes[1,1].plot(t_arr, results[obs]["hp"], label=f'obs z={x[obs[2]]:.1f}')
axes[1,1].set_xlabel('t'); axes[1,1].set_ylabel('h+')
axes[1,1].set_title('h+ at multiple observers')
axes[1,1].legend(); axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig12_quadrupole_gw.png", dpi=100)
print("Saved: paper/figures/fig12_quadrupole_gw.png")

np.savez("results/validated/v9_3_quadrupole_gw.npz",
         t=t_arr,
         hp_obs0=results[obs_list[0]]["hp"],
         hp_obs1=results[obs_list[1]]["hp"],
         hp_obs2=results[obs_list[2]]["hp"])
