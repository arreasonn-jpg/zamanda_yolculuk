
import numpy as np
from zt005.finite_wire_study import fields_with_groups
from zt005.model import Backpack,Drive

def test_group_sizes_and_current():
    pts=np.array([[.2,0,0]],float)
    _,_,_,sources,currents,groups=fields_with_groups("G1",pts,Backpack(),Drive(),.005,4,.01)
    assert len(groups)==6 and len(sources)==24
    for gi,g in enumerate(groups):
        st=sum(len(x) for x in groups[:gi])
        assert np.isclose(sum(currents[st:st+len(g)]),100.0)

def test_zero_radius_equals_one_filament_per_parent():
    pts=np.array([[.1,0,0]],float)
    _,_,T,sources,currents,groups=fields_with_groups("G1",pts,Backpack(),Drive(),0.,4,.01)
    assert len(groups)==6 and len(sources)==6
    assert np.allclose(currents,100.0)
    assert np.isfinite(T).all()
