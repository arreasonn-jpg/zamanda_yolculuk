
import numpy as np
from zt005.finite_wire_sources import thicken_source

def test_zero_radius_identity():
    pts=np.array([[1.,0.,0.],[1.,1.,0.]],float)
    seg=np.array([[0.,1.,0.],[0.,-1.,0.]],float)
    out=thicken_source((pts,seg),0.0,4)
    assert len(out)==1
    assert np.allclose(out[0][0],pts)

def test_positive_radius_expands_sources():
    pts=np.array([[1.,0.,0.],[1.,1.,0.]],float)
    seg=np.array([[0.,1.,0.],[0.,-1.,0.]],float)
    out=thicken_source((pts,seg),0.005,4)
    assert len(out)>=3
    assert any(not np.allclose(x[0],pts) for x in out)
