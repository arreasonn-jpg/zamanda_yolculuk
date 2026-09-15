
import numpy as np
from zt005.pointwise_curl import curl_point_2, curl_point_4

def test_pointwise_curl_constant_field():
    def field(p): return np.array([2.,-1.,3.])
    assert np.allclose(curl_point_2(field,[0.,0.,0.],0.01),0.)
    assert np.allclose(curl_point_4(field,[0.,0.,0.],0.01),0.)

def test_pointwise_curl_linear_field():
    # V=(-y,z,-x), curl(V)=(-1,1,1).
    def field(p):
        x,y,z=p
        return np.array([-y,z,-x],float)
    expected=np.array([-1.,1.,1.])
    assert np.allclose(curl_point_2(field,[0.,0.,0.],0.01),expected)
    assert np.allclose(curl_point_4(field,[0.,0.,0.],0.01),expected)
