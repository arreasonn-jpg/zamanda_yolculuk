"""v12.0-alpha2: Bowen-York binary BH orbit + GW extraction.
Two punctures with tangential momenta -> orbital motion.
Extract h+ at observer along z. Expected f_GW = 2 f_orbit.
"""
import numpy as np, time
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_ic_binary import gw_frequency
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.gw_extraction import extract_at_observer

N = 9; L = 20.0
m1 = 0.25; m2 = 0.25; d = 3.0
h = L/(N-1); dt = 2e-3; steps = 200

f_gw_expected = gw_frequency(m1, m2, d)
print(f"Binary orbit: N={N}, d={d}, f_GW_expected={f_gw_expected:.4f}")

s_base = bowen_york_binary(N=N, L=L, m1=m1, m2=m2, d=d)
s = BSSNStateMP.from_base(s_base)

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
mask = R > 2.0
H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
print(f"H0_ext = {H0:.4e}")

k_hat = np.array([0.0, 0.0, 1.0])
e1 = np.array([1.0, 0.0, 0.0])
e2 = np.array([0.0, 1.0, 0.0])
obs = (4, 4, 7)

t_hist = []; hp_hist = []; hc_hist = []; H_hist = []
t0 = time.time()
for step in range(1, steps+1):
    s = rk4_step_mp(s, h=h, dt=dt)
    if step % 5 == 0:
        hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
        t_hist.append(step*dt)
        hp_hist.append(hp); hc_hist.append(hc)
    if step % 50 == 0:
        H = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        H_hist.append(H)
        fin = (np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.K)))
        print(f"step={step:4d} H={H:.3e} growth={H/H0:.4f} "
              f"h+={hp_hist[-1]:.3e} amin={s.alpha.min():.4f} "
              f"fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break

t_arr = np.array(t_hist)
hp_arr = np.array(hp_hist); hc_arr = np.array(hc_hist)
print(f"\nh+ max = {np.max(np.abs(hp_arr)):.3e}")
print(f"hx max = {np.max(np.abs(hc_arr)):.3e}")
if len(H_hist) > 1:
    print(f"H growth = {H_hist[-1]/H0:.4f}")

if np.max(np.abs(hp_arr)) > 1e-20 and len(hp_arr) > 8:
    dt_out = t_arr[1] - t_arr[0]
    freqs = np.fft.rfftfreq(len(hp_arr), dt_out)
    spec = np.abs(np.fft.rfft(hp_arr))
    peak = freqs[1:][np.argmax(spec[1:])]
    print(f"FFT peak = {peak:.4f}, expected = {f_gw_expected:.4f}, "
          f"ratio = {peak/f_gw_expected:.3f}")

np.savez("results/validated/v12_0_binary_orbit.npz",
         t=t_arr, hp=hp_arr, hc=hc_arr, H=np.array(H_hist))
print("Saved results/validated/v12_0_binary_orbit.npz")
