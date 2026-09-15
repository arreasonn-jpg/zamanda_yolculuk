
import numpy as np
from zt005.stress_energy_conditioning import componentwise_errors, field_error_components

def test_component_errors_zero():
    T=np.ones((3,4,4)); E=np.ones((3,3)); B=np.ones((3,3))
    assert componentwise_errors(T,T)["T00"]["relative_rms"]==0
    assert field_error_components(E,B,E,B)["E"]["relative_l2"]==0
