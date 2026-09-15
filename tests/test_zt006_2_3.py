
import numpy as np
from zt005.error_budget import validation_points, relative_tensor_l2

def test_validation_points_inside():
    p=validation_points()
    assert p.ndim==2 and p.shape[1]==3
    assert np.all(p[:,0]**2+p[:,2]**2 <= .75**2+1e-12)
    assert np.all(np.abs(p[:,1]) <= 1.0+1e-12)

def test_relative_error_zero():
    a=np.ones((2,4,4))
    assert np.isclose(relative_tensor_l2(a,a),0.0)
