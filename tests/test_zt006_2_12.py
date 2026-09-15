
import numpy as np
from zt005.finite_wire_study import fields_with_groups
from zt005.model import Backpack,Drive

def test_parent_current_conservation():
    pts=np.array([[.2,0,0],[0,.2,0]],float)
    _,_,_,sources,currents,groups=fields_with_groups("G1",pts,Backpack(),Drive(),.005,4,.01)
    assert len(groups)==6 and len(sources)==24
    for gidx,g in enumerate(groups):
        start=sum(len(x) for x in groups[:gidx])
        assert np.isclose(sum(currents[start:start+len(g)]),Drive().peak_current_a)

def test_zero_radius_parent_current():
    pts=np.array([[.1,0,0]],float)
    _,_,T,sources,currents,groups=fields_with_groups("G1",pts,Backpack(),Drive(),0.,4,.01)
    assert len(groups)==6 and len(sources)==6
    assert np.allclose(currents,Drive().peak_current_a)
    assert np.isfinite(T).all()
