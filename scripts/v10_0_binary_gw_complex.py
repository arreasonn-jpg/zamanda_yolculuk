"""v10.0-alpha2: binary GW with complex scalar source.

Complex scalar (Psi_R + i Psi_I) with two peaks separated along
x-axis. T_ij has anisotropic content -> nonzero TT strain.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_complex import rk4_step_complex
from zt005.binary_complex_ic import complex_binary_ic
from zt005.gw_extraction import (
    extract_at_observer, psi4_from_history,
)

N = 11; L = 20.0; h = L/(N-1)
dt = 2e-3; steps = 200
A_R = 1e-3; A_I = 1e-3; sigma = 1.2; sep = 4.0; omega = 0.5

ic = complex_binary_ic(N=N, L=L, A_R=A_R, A_I=A_I, sigma=sigma,
                        separation=sep, omega=omega)
PR = ic["Psi_R"]; PI = ic["Psi_I"]
PiR = ic["Pi_R"]; PiI = ic["Pi_I"]
Z = np.zeros((N, N, N))

s = minkowski_state(N)
k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])
obs = (5, 5, 10)

print(f"Binary GW complex: omega={omega}, f_GW_expected={omega/np.pi:.4f}")

t_hist = []; hp_hist = []; hc_hist = []
for step in range(steps):
    s, PR, PI, PiR, PiI = rk4_step_complex(s, PR, PI, PiR, PiI, h, dt)
    t_hist.append((step+1)*dt)
    hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
    hp_hist.append(hp); hc_hist.append(hc)

t_arr = np.array(t_hist)
hp_arr = np.array(hp_hist); hc_arr = np.array(hc_hist)
dt_out = t_arr[1] - t_arr[0]

print(f"\nStrain: h+ max = {np.max(np.abs(hp_arr)):.3e}")
print(f"        hx max = {np.max(np.abs(hc_arr)):.3e}")

if np.max(np.abs(hp_arr)) > 1e-30:
    freqs = np.fft.rfftfreq(len(hp_arr), dt_out)
    spectrum = np.abs(np.fft.rfft(hp_arr))
    peak_idx = np.argmax(spectrum[1:]) + 1
    f_peak = freqs[peak_idx]
    print(f"FFT peak: {f_peak:.4f}")
    print(f"Expected f_GW = omega/pi = {omega/np.pi:.4f}")
    print(f"Ratio: {f_peak/(omega/np.pi):.3f}")
else:
    print("h+ is zero (no GW signal)")
    f_peak = 0; spectrum = np.zeros(1)

psi4 = psi4_from_history(hp_arr, hc_arr, dt_out)
print(f"psi4 max = {np.max(np.abs(psi4)) if len(psi4)>0 else 0:.3e}")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
axes[0].plot(t_arr, hp_arr, 'b-', label='h+')
axes[0].plot(t_arr, hc_arr, 'r--', label='hx')
axes[0].set_xlabel('t'); axes[0].set_ylabel('h')
axes[0].set_title('Binary GW strain (complex scalar)')
axes[0].legend(); axes[0].grid(True, alpha=0.3)

if np.max(np.abs(hp_arr)) > 1e-30:
    axes[1].plot(freqs[1:len(freqs)//2], spectrum[1:len(freqs)//2], 'r-')
    axes[1].axvline(omega/np.pi, color='k', ls='--',
                     label=f'expected f_GW={omega/np.pi:.3f}')
    axes[1].set_xlabel('f'); axes[1].set_ylabel('|FFT(h+)|')
    axes[1].set_title('GW frequency spectrum')
    axes[1].legend()

if len(psi4) > 0 and np.max(np.abs(psi4)) > 1e-30:
    t_psi = t_arr[:len(psi4)]
    axes[2].plot(t_psi, np.abs(psi4), 'g-')
    axes[2].set_xlabel('t'); axes[2].set_ylabel('|psi4|')
    axes[2].set_title('NP scalar psi4')

plt.tight_layout()
plt.savefig("paper/figures/fig13_binary_gw.png", dpi=100)
print("Saved: paper/figures/fig13_binary_gw.png")

np.savez("results/validated/v10_0_binary_gw.npz",
         t=t_arr, hp=hp_arr, hc=hc_arr,
         psi4=psi4 if len(psi4)>0 else np.array([]))
