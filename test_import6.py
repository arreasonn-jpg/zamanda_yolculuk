from zt005.bssn_ic_binary import brill_lindquist_state
from zt005.bssn_state_mp import BSSNStateMP
import numpy as np

s = brill_lindquist_state(33, 80.0, 0.5, 0.5, 10.0)
s_mp = BSSNStateMP.from_base(s)

# Check Atilde
A2 = np.einsum('...ij,...ij->...', s_mp.Atilde, s_mp.Atilde)
print('A2 shape:', A2.shape)
print('A2 max:', np.max(A2))
print('A2 min:', np.min(A2))
print('A2 mean:', np.mean(A2))

# Check profile along x
profile = A2.sum(axis=(1, 2))
print('Profile shape:', profile.shape)
print('Profile max:', np.max(profile))
print('Profile[:10]:', profile[:10])
print('Profile[-10:]:', profile[-10:])

# Check if there are any non-zero values
print('Non-zero count:', np.count_nonzero(A2))