
import numpy as np
from zt005.spatial_gr import inverse_trace_reverse, spatial_hbar_from_source

def test_inverse_trace_reverse_shape():
    hbar=np.zeros((4,4)); hbar[0,0]=2.
    h=inverse_trace_reverse(hbar)
    assert h.shape==(4,4)
    assert np.isfinite(h).all()

def test_spatial_green_linear_scaling():
    T=np.zeros((10,4,4)); T[:,0,0]=1.
    src=np.zeros((10,3)); obs=np.array([[1.,0.,0.]])
    h1=spatial_hbar_from_source(obs,src,T,1.0)
    h2=spatial_hbar_from_source(obs,src,T,2.0)
    assert np.allclose(h2,2*h1)

def test_no_nan_for_separated_points():
    T=np.zeros((5,4,4)); T[:,0,0]=1.
    src=np.array([[2.,0,0],[2.,1,0],[2.,0,1],[2.,-1,0],[2.,0,-1.]])
    obs=np.array([[0.,0,0]])
    h=spatial_hbar_from_source(obs,src,T,1.0)
    assert np.isfinite(h).all()
