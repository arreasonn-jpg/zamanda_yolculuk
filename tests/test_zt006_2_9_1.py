
import numpy as np
from zt005.gauss_cyl_quadrature import gauss_cylindrical_nodes

def test_gauss_volume_weights():
    pts,w=gauss_cylindrical_nodes(1.0,2.5,5,16,5)
    assert pts.shape==(5*16*5,3)
    assert w.shape==(5*16*5,)
    assert np.isclose(w.sum(),np.pi*2.5,rtol=1e-12)

def test_gauss_inside_cylinder():
    pts,_=gauss_cylindrical_nodes(1.0,2.5,4,8,4)
    assert np.all(pts[:,0]**2+pts[:,2]**2 <= 1.0+1e-12)
    assert np.all(np.abs(pts[:,1]) <= 1.25+1e-12)

def test_gauss_green_prefactor():
    from zt005.gauss_cyl_quadrature import integrate_green, EINSTEIN_GREEN_FACTOR
    pts=np.array([[1.,0.,0.]])
    vals=np.ones((1,4,4))
    w=np.array([1.])
    out=integrate_green(np.array([0.,0.,0.]),pts,vals,w)
    assert np.allclose(out,EINSTEIN_GREEN_FACTOR)
