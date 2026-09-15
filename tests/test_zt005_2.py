
from zt005.model import ActiveCylinder, Backpack
from zt005.production import evaluate_geometry

def test_geometry_fits():
    r=evaluate_geometry("G1",ActiveCylinder(),Backpack(),100.0,(17,17,17),64,96,2048)
    assert r["B_biot_savart_mean_T"] > 0
    assert r["B_curlA_mean_T"] > 0

def test_small_resolution_improves():
    a=evaluate_geometry("G1",ActiveCylinder(),Backpack(),100.0,(17,17,17),64,96,2048)
    b=evaluate_geometry("G1",ActiveCylinder(),Backpack(),100.0,(25,25,25),128,192,2048)
    assert b["B_curlA_vs_BS_relL2"] < a["B_curlA_vs_BS_relL2"]
