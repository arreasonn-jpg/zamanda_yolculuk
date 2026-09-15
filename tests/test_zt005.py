import numpy as np
from zt005.model import Backpack
from zt005.geometry import build_sources, check_source_envelope
from zt005.maxwell import vector_potential, biot_savart, curl_on_regular_grid

def test_all_sources_fit_backpack():
    b=Backpack()
    for g in ["G1","G2","G3","G4"]:
        assert check_source_envelope(build_sources(g,b),b)

def test_vector_potential_finite():
    b=Backpack()
    p=np.array([[0.0, b.center_y_m+0.2, 0.0]])
    A=vector_potential(p,build_sources("G1",b),100.0)
    assert A.shape==(1,3)
    assert np.all(np.isfinite(A))

def test_biot_savart_finite():
    b=Backpack()
    p=np.array([[0.0, b.center_y_m+0.2, 0.0]])
    B=biot_savart(p,build_sources("G1",b),100.0)
    assert np.all(np.isfinite(B))

def test_curl_shape():
    x=np.linspace(-1,1,5); y=np.linspace(-1,1,5); z=np.linspace(-1,1,5)
    X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    V=np.zeros(X.shape+(3,))
    V[...,2]=X
    C=curl_on_regular_grid(V,(x,y,z))
    assert C.shape==V.shape
    assert np.all(np.isfinite(C))
