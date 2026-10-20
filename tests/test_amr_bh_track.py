import numpy as np
from zt005.bssn_state_mp import BSSNStateMP
from zt005.bssn_ic_bowen_york import bowen_york_binary
from zt005.amr.amr_bh_track import find_bh_centers, pair_separation

def test_find_centers_initial():
    """Initial IC: punctures at x = +/- 2 (d=4)."""
    s = BSSNStateMP.from_base(bowen_york_binary(17, 20.0, 0.5, 0.5, 4.0))
    xl, xr, _ = find_bh_centers(s, L=20.0, N=17)
    sep = pair_separation(xl, xr)
    # Grid spacing 1.25, so 2 grid cells of error acceptable
    assert 2.0 <= sep <= 6.0, f"separation {sep:.2f} out of range"
