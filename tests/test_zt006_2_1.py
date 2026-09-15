
import numpy as np
from zt005.cyl_quadrature import cylindrical_midpoints, quadrature_integral

def test_cylindrical_volume():
    pts,w=cylindrical_midpoints(1.0,2.5,8,16,8)
    assert pts.shape[1]==3
    assert np.isclose(w.sum(),np.pi*1.0**2*2.5)

def test_constant_integrand():
    pts,w=cylindrical_midpoints(1.0,2.5,8,16,8)
    vals=np.ones((len(pts),4,4))
    out=quadrature_integral(np.array([0.,0.,0.]),pts,vals,w,softening_m=0.01)
    assert out.shape==(4,4)
    assert np.isfinite(out).all()

def test_einstein_green_prefactor():
    from zt005.cyl_quadrature import EINSTEIN_GREEN_FACTOR
    G=6.67430e-11
    C=299792458.0
    expected=4*G/C**4
    assert np.isclose(EINSTEIN_GREEN_FACTOR, expected)
