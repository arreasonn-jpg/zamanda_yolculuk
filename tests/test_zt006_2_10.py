
import numpy as np
from zt005.finite_wire_quadrature import min_distance_to_filament

def test_min_distance():
    pts=np.array([[0.,0.,0.],[2.,0.,0.]])
    src=[(np.array([[1.,0.,0.],[1.,0.,1.]]),np.zeros((2,3)))]
    d=min_distance_to_filament(pts,src)
    assert np.allclose(d,[1.,1.])

def test_finite_positive_distances():
    pts=np.array([[0.,0.,0.],[3.,0.,0.]])
    src=[(np.array([[1.,0.,0.],[1.,0.,1.]]),np.zeros((2,3)))]
    d=min_distance_to_filament(pts,src)
    assert np.all(d>0)
