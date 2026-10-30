with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix the import - use sommerfeld_all instead of sommerfeld_boundary
old_import = """from zt005.bssn_ic_binary import brill_lindquist_state as brill_lindquist_ic
from zt005.bssn_state import BSSNState
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
from zt005.bssn_dissipation import sommerfeld_boundary as apply_sommerfeld_boundary
from zt005.bssn_rk4_beta3 import rk4_step_beta3"""

new_import = """from zt005.bssn_ic_binary import brill_lindquist_state as brill_lindquist_ic
from zt005.bssn_state import BSSNState
from zt005.bssn_constraints import hamiltonian_constraint_grid
from zt005.amr.amr_projection import project_phi_band_covariant, momentum_band_correct, _band_mask
from zt005.bssn_bc import sommerfeld_all
from zt005.bssn_rk4_beta3 import rk4_step_beta3"""

content = content.replace(old_import, new_import)

# Replace the apply_sommerfeld_boundary calls with sommerfeld_all
# Need to keep track of old state for sommerfeld_all
content = content.replace(
    'state_c = rk4_step_beta3(state_c, dx_c, dt)\n        apply_sommerfeld_boundary(state_c)',
    'old_c = BSSNState(\n            phi=state_c.phi.copy(), gtilde=state_c.gtilde.copy(),\n            K=state_c.K.copy(), Atilde=state_c.Atilde.copy(),\n            Gamma=state_c.Gamma.copy(), alpha=state_c.alpha.copy(),\n            beta=state_c.beta.copy()\n        )\n        state_c = rk4_step_beta3(state_c, dx_c, dt)\n        bc_c = sommerfeld_all(state_c, old_c, dx_c, dt)\n        state_c = BSSNState(\n            phi=bc_c["phi"], gtilde=bc_c["gtilde"],\n            K=bc_c["K"], Atilde=bc_c["Atilde"],\n            Gamma=bc_c["Gamma"], alpha=bc_c["alpha"],\n            beta=bc_c["beta"]\n        )'
)

content = content.replace(
    'state_f = rk4_step_beta3(state_f, dx_f, dt)\n        apply_sommerfeld_boundary(state_f)',
    'old_f = BSSNState(\n            phi=state_f.phi.copy(), gtilde=state_f.gtilde.copy(),\n            K=state_f.K.copy(), Atilde=state_f.Atilde.copy(),\n            Gamma=state_f.Gamma.copy(), alpha=state_f.alpha.copy(),\n            beta=state_f.beta.copy()\n        )\n        state_f = rk4_step_beta3(state_f, dx_f, dt)\n        bc_f = sommerfeld_all(state_f, old_f, dx_f, dt)\n        state_f = BSSNState(\n            phi=bc_f["phi"], gtilde=bc_f["gtilde"],\n            K=bc_f["K"], Atilde=bc_f["Atilde"],\n            Gamma=bc_f["Gamma"], alpha=bc_f["alpha"],\n            beta=bc_f["beta"]\n        )'
)

with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('Fixed Sommerfeld boundary usage')