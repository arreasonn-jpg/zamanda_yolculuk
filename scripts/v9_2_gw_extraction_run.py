"""v9.2-alpha2: BSSN evolution with GW extraction.

NULL RESULT: Spherically symmetric scalar perturbation produces
zero TT strain (h+ = hx = 0). This is physically correct:
a monopole source has no quadrupole moment, so no gravitational
waves are emitted. The BSSN + NP pipeline works, but this
specific setup has no GW signal.

To produce GW, a quadrupole source is required (future work).


Setup: Minkowski + small scalar perturbation
    phi_sc ~ 1e-4 Gaussian

Evolve 200 steps with bssn_rk4_matter.
Extract h_ij at 3 observer positions.
Compute psi_4 from time history.
Recover |h| strain amplitude.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.gw_extraction import (
    metric_perturbation_from_state, extract_at_observer,
    psi4_from_history, strain_amplitude_from_psi4_history,
)

N = 9; L = 20.0; h = L/(N-1); dt = 1e-3; steps = 50
A_pert = 1e-4; sigma = 1.5

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R2 = X**2 + Y**2 + Z**2

# Initial scalar Gaussian
phi_sc = A_pert * np.exp(-R2 / (2*sigma**2))
Pi_sc = np.zeros((N, N, N))
V_func = lambda p: 0.5 * p**2

s = minkowski_state(N)

# Observer direction (z)
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])

# Observer indices
obs_list = [(4, 4, 6), (4, 4, 7), (4, 4, 8)]  # along z

hp_hist = {obs: [] for obs in obs_list}
hc_hist = {obs: [] for obs in obs_list}
t_hist = []

print(f"BSSN evolution with GW extraction: N={N}, steps={steps}")

for step in range(steps):
    s, phi_sc, Pi_sc = rk4_step_matter(s, phi_sc, Pi_sc, h, dt, V_func=V_func)
    if step % 2 == 0:
        t_hist.append((step+1)*dt)
        for obs in obs_list:
            hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
            hp_hist[obs].append(hp)
            hc_hist[obs].append(hc)

t_arr = np.array(t_hist)
print(f"Time series: {len(t_arr)} samples, dt_out = {t_arr[1]-t_arr[0]:.2e}")

# Compute psi_4 for each observer
results = {}
dt_out = t_arr[1] - t_arr[0]
for obs in obs_list:
    hp = np.array(hp_hist[obs]); hc = np.array(hc_hist[obs])
    psi4 = psi4_from_history(hp, hc, dt_out)
    h_est = strain_amplitude_from_psi4_history(psi4, dt_out)
    results[obs] = {"hp": hp, "hc": hc, "psi4": psi4, "h_est": h_est}
    print(f"obs={obs}: h+ max={np.max(np.abs(hp)):.3e}, "
          f"psi4 max={np.max(np.abs(psi4)) if len(psi4)>0 else 0:.3e}, "
          f"h_est max={np.max(np.abs(h_est)) if len(h_est)>0 else 0:.3e}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

obs0 = obs_list[-1]  # farthest observer
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
    hp = results[obs]["hp"]
    axes[1,1].plot(t_arr, hp, label=f'obs z={x[obs[2]]:.1f}')
axes[1,1].set_xlabel('t'); axes[1,1].set_ylabel('h+')
axes[1,1].set_title('h+ at multiple observers')
axes[1,1].legend(); axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig11_gw_extraction.png", dpi=100)
print("Saved: paper/figures/fig11_gw_extraction.png")

np.savez("results/validated/v9_2_gw_extraction.npz",
         t=t_arr,
         hp_obs0=results[obs_list[0]]["hp"],
         hp_obs1=results[obs_list[1]]["hp"],
         hp_obs2=results[obs_list[2]]["hp"])
