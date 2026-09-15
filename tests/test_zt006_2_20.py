import numpy as np
from zt005.external_green_validation import external_observers
def test_observers_outside():
    assert np.all(np.linalg.norm(external_observers(),axis=1)>1.0)
def test_count():
    assert len(external_observers())==5
