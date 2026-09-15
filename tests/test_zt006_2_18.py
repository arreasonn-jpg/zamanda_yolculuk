
import numpy as np
from zt005.singular_green import regularized_green_integral

def test_regularized_integral_finite():
    p=np.array([[0.,0.,0.],[1.,0.,0.]])
    w=np.array([1.,1.])
    T=np.zeros((2,4,4)); T[:,0,0]=1
    h=regularized_green_integral(np.array([0.,0.,0.]),p,w,T,1.0,0.5)
    assert np.isfinite(h).all()
    assert h[0,0] > 0

def test_zero_near_is_finite():
    p=np.array([[0.,0.,0.]])
    w=np.array([1.])
    T=np.ones((1,4,4))
    h=regularized_green_integral(np.array([0.,0.,0.]),p,w,T,1.0,0.1)
    assert np.isfinite(h).all()
