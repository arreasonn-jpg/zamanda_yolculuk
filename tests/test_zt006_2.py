import numpy as np
from zt005.scaling_study import make_source_points
def test_source_points_inside():
    class C: radius_m=1.0; height_m=2.5
    p=make_source_points(C(),200)
    assert p.shape==(200,3)
    assert np.all(p[:,0]**2+p[:,2]**2 <= .95**2+1e-12)
    assert np.all(np.abs(p[:,1]) <= .95*2.5/2 + 1e-12)
