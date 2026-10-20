import zt005.amr.amr_interp
import zt005.amr.amr_bh_track
import zt005.amr.refinement
import zt005.amr.refinement_tensor
import zt005.amr.bssn_amr_consistent
import zt005.matter.matter_scalar
import zt005.matter.matter_complex
import zt005.matter.matter_ptum
import zt005.matter.matter_proca
import zt005.matter.matter_brans_dicke
import zt005.matter.sc_condensate
import zt005.matter.exotic_matter
import zt005.matter.phase_transition_common
import zt005.matter.ptum_field

# old-style paths must also still resolve (repr via sys.modules keys)
import sys
mods = [name for name in sys.modules if name.startswith("zt005.")]
print("smoke ok, modules:", len(mods))
