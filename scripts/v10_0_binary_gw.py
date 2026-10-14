"""v10.0-alpha1: binary orbit GW extraction.

Two scalar Gaussian peaks in circular orbit around COM.
Extract h+ at observer; expect f_GW = 2 * f_orbit signal.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter
from zt005.binary_orbit_ic import (
    gaussian_binary_ic, expected_gw_frequency,
)
from zt005.gw_extraction import (
    extract_at_observer, psi4_from_history,
)

N = 11; L = 20.0; h = L/(N-1)
dt = 2e-3; steps = 300
A_amp = 1e-3; sigma = 1.2; sep = 4.0
omega = 0.5  # rad / time

ic = gaussian_binary_ic(N=N, L=L, A_amp=A_amp, sigma=sigma,
                         separation=sep, omega=omega)
phi_sc = ic["phi_sc"]; Pi_sc = ic["Pi_sc"]
V_func = lambda p: 0.5 * p**2

s = minkowski_state(N)
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])
obs = (5, 5, 10)

print(f"Binary GW: N={N}, omega={omega}, f_orbit={omega/(2*np.pi):.4f}, "
      f"f_GW_expected={expected_gw_frequency(omega):.4f}")

t_hist = []; hp_hist = []; hc_hist = []
for step in range(steps):
    s, phi_sc, Pi_sc = rk4_step_matter(s, phi_sc, Pi_sc, h, dt, V_func=V_func)
    t_hist.append((step+1)*dt)
    hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
    hp_hist.append(hp); hc_hist.append(hc)

t_arr = np.array(t_hist)
hp_arr = np.array(hp_hist); hc_arr = np.array(hc_hist)
dt_out = t_arr[1] - t_arr[0]

# FFT for frequency
freqs = np.fft.rfftfreq(len(hp_arr), dt_out)
spectrum = np.abs(np.fft.rfft(hp_arr))
peak_idx = np.argmax(spectrum[1:]) + 1
f_peak = freqs[peak_idx]

print(f"\nStrain: h+ max = {np.max(np.abs(hp_arr)):.3e}")
print(f"FFT peak frequency: {f_peak:.4f}")
print(f"Expected f_GW = 2*f_orbit = {expected_gw_frequency(omega):.4f}")
print(f"Ratio: f_peak / f_GW_expected = {f_peak/expected_gw_frequency(omega):.3f}")

psi4 = psi4_from_history(hp_arr, hc_arr, dt_out)
print(f"psi4 max = {np.max(np.abs(psi4)) if len(psi4)>0 else 0:.3e}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
axes[0].plot(t_arr, hp_arr, 'b-')
axes[0].set_xlabel('t'); axes[0].set_ylabel('h+')
axes[0].set_title('Binary orbit GW strain')
axes[0].grid(True, alpha=0.3)

axes[1].plot(freqs[1:len(freqs)//2], spectrum[1:len(freqs)//2], 'r-')
axes[1].axvline(expected_gw_frequency(omega), color='k', ls='--',
                 label=f'expected f_GW={expected_gw_frequency(omega):.3f}')
axes[1].set_xlabel('f'); axes[1].set_ylabel('|FFT(h+)|')
axes[1].set_title('GW frequency spectrum')
axes[1].legend(); axes[1].grid(True, alpha=0.3)

if len(psi4) > 0:
    t_psi = t_arr[:len(psi4)]
    axes[2].plot(t_psi, np.abs(psi4), 'g-')
    axes[2].set_xlabel('t'); axes[2].set_ylabel('|psi4|')
    axes[2].set_title('NP scalar psi4')
    axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("paper/figures/fig13_binary_gw.png", dpi=100)
print("Saved: paper/figures/fig13_binary_gw.png")

np.savez("results/validated/v10_0_binary_gw.npz",
         t=t_arr, hp=hp_arr, hc=hc_arr,
         freqs=freqs, spectrum=spectrum, f_peak=f_peak)
