from zt005.bssn_ic_binary import brill_lindquist_state
import numpy as np
s = brill_lindquist_state(33, 80.0, 0.5, 0.5, 10.0)
print('OK:', s.phi.shape)