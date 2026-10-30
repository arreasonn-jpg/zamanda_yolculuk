from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_rk4_mp_vec_damp import rk4_step_mp_vec_damp
from zt005.amr.refinement_tensor import prolongate_state
from zt005.amr.amr_interp import inject_boundary
from zt005.amr.refinement import prolongate
from zt005.amr.amr_moving import fine_grid_mask, target_lo
from zt005.amr.amr_bh_track import peak_centroids, pair_separation
import numpy as np

# Test update_fine_patch
s = bowen_york_binary(33, 80.0, 0.5, 0.5, 10.0)
s_c = BSSNStateMP.from_base(s)
print('Coarse OK:', s_c.phi.shape)

# Test peak_centroids
c_left, c_right, _ = peak_centroids(s_c, 80.0)
print('Peak centroids:', c_left, c_right)
com_x = (c_left[0] + c_right[0]) / 2.0
print('COM:', com_x)

# Test target_lo
h_c = 80.0 / 32
lo = target_lo(com_x, h_c, 33, 25)
print('lo:', lo)

# Test update_fine_patch logic
new_offset = (lo, lo, lo)
fine_dict2 = prolongate_state(s_c, 25, new_offset)
B_f2 = np.zeros((25, 25, 25, 3))
for i in range(3):
    B_f2[..., i] = prolongate(s_c.B[..., i], 25, new_offset)
fine_dict2["B"] = B_f2

coarse_prolongated = prolongate_state(s_c, 25, new_offset)
for key in ["phi", "K", "alpha", "gtilde", "Atilde", "Gamma", "beta"]:
    fine_dict2[key] = inject_boundary(
        fine_dict2[key], coarse_prolongated[key], n_layers=4
    )

s_f_new = BSSNStateMP(
    phi=fine_dict2["phi"],
    gtilde=fine_dict2["gtilde"],
    K=fine_dict2["K"],
    Atilde=fine_dict2["Atilde"],
    Gamma=fine_dict2["Gamma"],
    alpha=fine_dict2["alpha"],
    beta=fine_dict2["beta"],
    B=fine_dict2["B"],
)
print('New fine OK:', s_f_new.phi.shape)

# Test RK4 on fine
h_f = 40.0 / 24
s_f_new2 = rk4_step_mp_vec_damp(s_f_new, h_f, 0.03125, eta=1.0, width_cells=3, kappa=2.0)
print('Fine RK4 OK:', s_f_new2.phi.shape)

print('All tests passed!')