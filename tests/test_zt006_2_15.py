
import numpy as np
from zt005.area_weighted_finite_wire import gauss_disk_nodes
def test_weights():
    p,w=gauss_disk_nodes(3,8)
    assert p.shape==(24,2)
    assert np.isclose(w.sum(),1.0)
def test_inside():
    p,_=gauss_disk_nodes(4,10)
    assert np.all(np.linalg.norm(p,axis=1)<=1+1e-12)
