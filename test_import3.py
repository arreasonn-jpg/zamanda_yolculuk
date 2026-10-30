from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_state_mp import BSSNStateMP
from zt005.amr.refinement_tensor import prolongate_state
import numpy as np
s = brill_lindquist_state(33, 80.0, 0.5, 0.5, 10.0)
s_mp = BSSNStateMP.from_base(s)
print('MP OK:', s_mp.phi.shape)
offset = ((33 - 25) // 2,) * 3
print('offset:', offset)
fine_dict = prolongate_state(s_mp, 25, offset)
print('Prolongation OK:', fine_dict['phi'].shape)