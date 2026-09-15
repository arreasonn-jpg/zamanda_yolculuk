
import numpy as np
from zt005.validated_stress_energy import em_stress_energy

def test_zero_field_zero_tensor():
    E=np.zeros((2,3)); B=np.zeros((2,3))
    s=em_stress_energy(E,B)
    assert np.allclose(s["u"],0)
    assert np.allclose(s["S"],0)
    assert np.allclose(s["T"],0)

def test_positive_energy_density():
    E=np.array([[1.,0.,0.]])
    B=np.array([[0.,1e-6,0.]])
    s=em_stress_energy(E,B)
    assert s["u"][0] > 0
    assert np.isfinite(s["T"]).all()
