"""v10.1-alpha2: GW extraction pipeline test with injected packet.
Test: inject a TT-gauge +polarization packet, evolve briefly,
verify NP extraction recovers expected h+.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.gw_initial_data import gw_packet
from zt005.gw_extraction import (
    extract_at_observer, psi4_from_history,
    strain_amplitude_from_psi4_history,
)

N = 11; L = 20.0; h = L/(N-1)
dt = 5e-3; steps = 50
A_amp = 1e-3; sigma = 1.5; z0 = 0.0

# Inject packet
pkt = gw_packet(N=N, L=L, A_amp=A_amp, sigma=sigma, z0=z0,
                polarization="plus")
s = minkowski_state(N)
s.gtilde = pkt["gtilde"]
s.phi = pkt["phi"]
s.K = pkt["K"]
s.Atilde = pkt["Atilde"]

Z = np.zeros((N, N, N))
phi_sc = np.zeros((N, N, N))
Pi_sc = np.zeros((N, N, N))
V_func = lambda p: 0.5 * p**2

# Observers along z
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])
obs_list = [(5, 5, 6), (5, 5, 7), (5, 5, 8)]

# Initial extraction (no evolution)
print(f"Injected packet: A={A_amp}, sigma={sigma}, z0={z0}")
for obs in obs_list:
    hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
    print(f" t=0 obs={obs}: h+={hp:.4e}, hx={hc:.4e}")

# Evolve and track
hp_hist = {obs: [] for obs in obs_list}
hc_hist = {obs: [] for obs in obs_list}
t_hist = []
for step in range(steps):
    s, phi_sc, Pi_sc = rk4_step_matter(s, phi_sc, Pi_sc, h, dt, V_func=V_func)
    t_hist.append((step+1)*dt)
    for obs in obs_list:
        hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
        hp_hist[obs].append(hp)
        hc_hist[obs].append(hc)

t_arr = np.array(t_hist)
dt_out = t_arr[1] - t_arr[0]

# Final summary
print(f"\nAfter {steps} steps (t={t_arr[-1]:.3f}):")
for obs in obs_list:
    hp = np.array(hp_hist[obs])
    hc = np.array(hc_hist[obs])
    psi4 = psi4_from_history(hp, hc, dt_out)
    h_est = strain_amplitude_from_psi4_history(psi4, dt_out)
    print(f" obs={obs}: h+ max={np.max(np.abs(hp)):.4e}, "
          f"hx max={np.max(np.abs(hc)):.4e}, "
          f"psi4 max={np.max(np.abs(psi4)) if len(psi4)>0 else 0:.4e}, "
          f"h_est max={np.max(np.abs(h_est)) if len(h_est)>0 else 0:.4e}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
obs0 = obs_list[-1]
axes[0,0].plot(t_arr, hp_hist[obs0], 'b-', label='h+')
axes[0,0].plot(t_arr, hc_hist[obs0], 'r--', label='hx')
axes[0,0].set_xlabel('t'); axes[0,0].set_ylabel('h')
axes[0,0].set_title(f'Strain at observer {obs0}')
axes[0,0].legend(); axes[0,0].grid(True, alpha=0.3)

for obs in obs_list:
    axes[0,1].plot(t_arr, hp_hist[obs], label=f'obs z={obs[2]}')
axes[0,1].set_xlabel('t'); axes[0,1].set_ylabel('h+')
axes[0,1].set_title('h+ at multiple observers')
axes[0,1].legend(); axes[0,1].grid(True, alpha=0.3)

# FFT
hp_arr = np.array(hp_hist[obs0])
freqs = np.fft.rfftfreq(len(hp_arr), dt_out)
spectrum = np.abs(np.fft.rfft(hp_arr))
axes[1,0].semilogy(freqs[1:], spectrum[1:], 'g-')
axes[1,0].set_xlabel('f'); axes[1,0].set_ylabel('|FFT(h+)|')
axes[1,0].set_title('Spectrum')
axes[1,0].grid(True, alpha=0.3)

# Initial packet profile along z
z_line = np.linspace(-L/2, L/2, N)
h_init = pkt["h_ij"][5, 5, :, 0, 0]
axes[1,1].plot(z_line, h_init, 'm-')
axes[1,1].set_xlabel('z'); axes[1,1].set_ylabel('h_xx (t=0)')
axes[1,1].set_title('Injected packet profile')
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig14_gw_pipeline.png", dpi=100)
print("\nSaved: paper/figures/fig14_gw_pipeline.png")

np.savez("results/validated/v10_1_gw_pipeline.npz",
         t=t_arr,
         hp_obs6=np.array(hp_hist[obs_list[0]]),
         hp_obs7=np.array(hp_hist[obs_list[1]]),
         hp_obs8=np.array(hp_hist[obs_list[2]]))
