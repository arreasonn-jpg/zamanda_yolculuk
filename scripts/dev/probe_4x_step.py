"""Profile a single 4x-AMR step (find the hot path, no artifacts)."""
import cProfile
import pstats
import io
import time

import numpy as np

from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_rk4_mp_vec_z4c import rk4_step_mp_vec_z4c
from zt005.amr.amr_interp import (
    prolongate_scalar, prolongate_tensor, prolongate_vector)

N_c, L_c = 17, 20.0
N_sub = 9
factor = 4
N_f = (N_sub - 1) * factor + 1
h_c = L_c / (N_c - 1)
L_f = 10.0
h_f = L_f / (N_f - 1)

s_c = BSSNStateMP.from_base(bowen_york_binary(N_c, L_c, 0.5, 0.5, 4.0))
lo, hi = 4, 13
win = slice(lo, hi)


def prolong(state):
    return BSSNStateMP(
        phi=prolongate_scalar(s_c.phi[win, win, win], factor),
        gtilde=prolongate_tensor(s_c.gtilde[win, win, win], factor),
        K=prolongate_scalar(s_c.K[win, win, win], factor),
        Atilde=prolongate_tensor(s_c.Atilde[win, win, win], factor),
        Gamma=prolongate_vector(s_c.Gamma[win, win, win], factor),
        alpha=prolongate_scalar(s_c.alpha[win, win, win], factor),
        beta=prolongate_vector(s_c.beta[win, win, win], factor),
        B=prolongate_vector(s_c.B[win, win, win], factor),
    )


t0 = time.time()
s_f = prolong(s_c)
print(f"prolong time: {time.time()-t0:.3f}s  (N_f={N_f})")

pr = cProfile.Profile()
pr.enable()
t0 = time.time()
for _ in range(3):
    s_c2 = rk4_step_mp_vec_z4c(s_c, h=h_c, dt=2e-3, width_cells=3,
                               kappa=0.0, eta=1.0)
t_c = time.time() - t0

t0 = time.time()
for _ in range(3):
    s_f2 = rk4_step_mp_vec_z4c(s_f, h=h_f, dt=5e-4, width_cells=3,
                               kappa=0.0, eta=1.0)
t_f = time.time() - t0
pr.disable()
print(f"coarse RK4 step: {t_c/3:.4f}s   fine RK4 step (N=33): {t_f/3:.4f}s"
      f"   -> x{t_f / max(t_c, 1e-9):.1f}")

s = io.StringIO()
stats = pstats.Stats(pr, stream=s).sort_stats("cumulative")
stats.print_stats(18)
print("\n".join(s.getvalue().splitlines()[:30]))
