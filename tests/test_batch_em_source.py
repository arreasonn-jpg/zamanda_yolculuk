
import numpy as np
from zt005.batch_em_source import compute_T_batch
from zt005.model import Backpack, Drive

def test_batch_tensor_shape():
    points=np.array([[0.2,0.0,0.0],[-0.2,0.1,0.1]])
    T=compute_T_batch("G1",points,Backpack(),Drive(),h=0.01)
    assert T.shape==(2,4,4)
    assert np.isfinite(T).all()
