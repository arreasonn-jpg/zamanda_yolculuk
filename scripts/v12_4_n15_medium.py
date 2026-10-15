"""v12.4-alpha1: N=15 medium run (1000 steps)."""
import numpy as np, time
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_ic_binary import gw_frequency
from zt005.bssn_rk4_mp_vec import rk4_step_mp_vec
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.gw_extraction import extract_at_observer

N = 15; L = 20.0
m1 = m2 = 0.5; d = 4.0
h = L/(N-1); dt = 1e-3; steps = 1000
f_gw = gw_frequency(m1, m2, d)
print(f"N={N} f_GW_expected={f_gw:.4f} steps={steps}")

s = BSSNStateMP.from_base(bowen_york_binary(N, L, m1, m2, d))
x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
mask = np.sqrt(X**2+Y**2+Z**2) > 1.5
H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))

k_hat = np.array([0., 0., 1.]); e1 = np.array([1., 0., 0.]); e2 = np.array([0., 1., 0.])
obs = (7, 7, 12)
t_hist, hp_hist, hc_hist = [], [], []
t0 = time.time()
for step in range(1, steps+1):
    s = rk4_step_mp_vec(s, h=h, dt=dt)
    if step % 5 == 0:
        hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
        t_hist.append(step*dt); hp_hist.append(hp); hc_hist.append(hc)
    if step % 200 == 0:
        H = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        fin = np.all(np.isfinite(s.phi)) and np.all(np.isfinite(s.K))
        print(f"step={step:4d} H/H0={H/H0:.4f} h+={hp_hist[-1]:.3e} "
              f"hx={hc_hist[-1]:.3e} fin={fin} t={time.time()-t0:.0f}s")
        if not fin:
            print("NON-FINITE"); break
t_arr = np.array(t_hist); hp = np.array(hp_hist); hc = np.array(hc_hist)
print(f"\nh+ max={np.max(np.abs(hp)):.3e} hx max={np.max(np.abs(hc)):.3e}")
np.savez("results/validated/v12_4_n15_medium.npz", t=t_arr, hp=hp, hc=hc)
