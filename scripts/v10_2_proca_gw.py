"""v10.2-alpha3: Proca vector field GW extraction.
Setup: E1-polarized vector wave packet propagating along z.
Expectation: anisotropic T_ij -> nonzero h+ at observer
(unlike scalar case in v10.0).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_proca import rk4_step_proca
from zt005.gw_extraction import extract_at_observer

N = 11; L = 20.0; h = L/(N-1)
dt = 5e-3; steps = 100
sigma = 1.5; z0 = 5.0; A_amp = 1e-3

x = np.linspace(-L/2, L/2, N)
X, Y, Z2 = np.meshgrid(x, x, x, indexing="ij")
env = A_amp * np.exp(-((Z2 - z0)**2) / (2.0 * sigma**2))

# Proca wave packet: E1 nonzero
E1 = env.copy()
A1 = np.zeros_like(E1)
Z = np.zeros((N, N, N))

s = minkowski_state(N)
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])
obs_list = [(5, 5, 6), (5, 5, 7), (5, 5, 8)]

hp_hist = {o: [] for o in obs_list}
t_hist = []
for step in range(steps):
    s, A1, Z, Z, E1, Z, Z = rk4_step_proca(
        s, A1, Z, Z, E1, Z, Z, h, dt, m2=0.0)
    t_hist.append((step + 1) * dt)
    for o in obs_list:
        hp, hc = extract_at_observer(s, o, k_hat, e1, e2)
        hp_hist[o].append(hp)

t_arr = np.array(t_hist)
print(f"Proca GW (E1 packet, z0={z0}):")
for o in obs_list:
    hp = np.array(hp_hist[o])
    print(f" obs={o}: h+ max={np.max(np.abs(hp)):.4e}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for o in obs_list:
    axes[0].plot(t_arr, hp_hist[o], label=f'obs z={o[2]}')
axes[0].set_xlabel('t'); axes[0].set_ylabel('h+')
axes[0].set_title('Proca GW strain (E1 packet)')
axes[0].legend(); axes[0].grid(True, alpha=0.3)

axes[1].plot(x, np.abs(E1[5, 5, :]), 'r-', label='|E1| final')
axes[1].plot(x, np.abs(A1[5, 5, :]), 'b-', label='|A1| final')
axes[1].set_xlabel('z'); axes[1].set_ylabel('amplitude')
axes[1].set_title('Proca field profile')
axes[1].legend(); axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig15_proca_gw.png", dpi=100)
print("Saved: paper/figures/fig15_proca_gw.png")
np.savez("results/validated/v10_2_proca_gw.npz",
         t=t_arr,
         hp_obs6=np.array(hp_hist[obs_list[0]]),
         hp_obs7=np.array(hp_hist[obs_list[1]]),
         hp_obs8=np.array(hp_hist[obs_list[2]]))
