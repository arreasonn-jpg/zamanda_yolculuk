"""v12.1-alpha1: Brill-Lindquist plunge GW waveform.
Longer run to capture the infall -> merger-like burst.
N=11, 400 steps.
"""
import numpy as np
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_rk4_mp import rk4_step_mp
from zt005.bssn_ricci_full import hamiltonian_constraint_full
from zt005.gw_extraction import extract_at_observer

N = 7; L = 20.0
m1 = 0.5; m2 = 0.5; d = 4.0
h = L/(N-1); dt = 5e-3; steps = 50

s = BSSNStateMP.from_base(brill_lindquist_state(N, L, m1, m2, d))

x = np.linspace(-L/2, L/2, N)
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
mask = np.sqrt(X**2+Y**2+Z**2) > 2.0
H0 = float(np.mean(np.abs(hamiltonian_constraint_full(
    s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))

k_hat = np.array([0., 0., 1.]); e1 = np.array([1., 0., 0.]); e2 = np.array([0., 1., 0.])
obs = (3, 3, 5)
t_hist, hp_hist, hc_hist = [], [], []
for step in range(1, steps+1):
    s = rk4_step_mp(s, h=h, dt=dt)
    if step % 2 == 0:
        hp, hc = extract_at_observer(s, obs, k_hat, e1, e2)
        t_hist.append(step*dt); hp_hist.append(hp); hc_hist.append(hc)
    if step % 50 == 0:
        H = float(np.mean(np.abs(hamiltonian_constraint_full(
            s.phi, s.gtilde, s.K, s.Atilde, h)[mask])))
        print(f"step={step:4d} H/H0={H/H0:.4f} h+={hp_hist[-1]:.3e}")
t_arr = np.array(t_hist); hp = np.array(hp_hist); hc = np.array(hc_hist)
print(f"h+ max = {np.max(np.abs(hp)):.3e}, hx max = {np.max(np.abs(hc)):.3e}")
np.savez("results/validated/v12_1_bl_plunge.npz", t=t_arr, hp=hp, hc=hc)
