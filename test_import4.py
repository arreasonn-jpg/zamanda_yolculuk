from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_rk4_mp_vec_damp import rk4_step_mp_vec_damp
import numpy as np
s = brill_lindquist_state(33, 80.0, 0.5, 0.5, 10.0)
s_mp = BSSNStateMP.from_base(s)
print('MP OK:', s_mp.phi.shape)
h = 80.0 / 32
print('h:', h)
s_new = rk4_step_mp_vec_damp(s_mp, h, 0.03125, eta=1.0, width_cells=3, kappa=2.0)
print('RK4 OK:', s_new.phi.shape)