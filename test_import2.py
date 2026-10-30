from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_state_mp import BSSNStateMP
import numpy as np
s = brill_lindquist_state(33, 80.0, 0.5, 0.5, 10.0)
print('OK:', s.phi.shape)
s_mp = BSSNStateMP.from_base(s)
print('MP OK:', s_mp.phi.shape)