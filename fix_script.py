with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace the custom rk4_step with proper import and usage
old_import = """from zt005.bssn_ic_binary import brill_lindquist_state as brill_lindquist_ic
from zt005.bssn_rhs_vec import bssn_rhs_vec
from zt005.bssn_state import BSSNState
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
from zt005.bssn_dissipation import sommerfeld_boundary as apply_sommerfeld_boundary

def rk4_step(state, rhs_fn, dt, *args, **kwargs):
    \"\"\"RK4 integrator.\"\"\"
    k1 = rhs_fn(state, *args, **kwargs)
    k2 = rhs_fn(state + 0.5 * dt * k1, *args, **kwargs)
    k3 = rhs_fn(state + 0.5 * dt * k2, *args, **kwargs)
    k4 = rhs_fn(state + dt * k3, *args, **kwargs)
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)"""

new_import = """from zt005.bssn_ic_binary import brill_lindquist_state as brill_lindquist_ic
from zt005.bssn_state import BSSNState
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
from zt005.bssn_dissipation import sommerfeld_boundary as apply_sommerfeld_boundary
from zt005.bssn_rk4_beta3 import rk4_step_beta3"""

content = content.replace(old_import, new_import)

# Replace the rk4_step calls with rk4_step_beta3
content = content.replace(
    'state_c = rk4_step(state_c, lambda s: bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, dx_c, kappa=kappa, beta=s.beta, Gamma=s.Gamma), dt)',
    'state_c = rk4_step_beta3(state_c, dx_c, dt)'
)
content = content.replace(
    'state_f = rk4_step(state_f, lambda s: bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, dx_f, kappa=kappa, beta=s.beta, Gamma=s.Gamma), dt)',
    'state_f = rk4_step_beta3(state_f, dx_f, dt)'
)

with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('Fixed RK4 step usage')