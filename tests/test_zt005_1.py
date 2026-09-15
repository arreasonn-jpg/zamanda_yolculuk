
import numpy as np
from zt005.convergence import analytic_loop_axis_B, benchmark_single_loop

def test_analytic_loop_axis_formula_positive_and_symmetric():
    b=analytic_loop_axis_B(np.array([-0.2,0.0,0.2]),0.1,100.0)
    assert b[1] < 0
    assert np.isclose(b[0],b[2])

def test_filament_convergence_improves():
    d=benchmark_single_loop(segment_counts=(32,128,512),grid_counts=(17,))
    e=[x["relative_L2_error"] for x in d["segment_convergence"]]
    assert e[-1] < e[0]

def test_grid_convergence_improves():
    d=benchmark_single_loop(segment_counts=(128,),grid_counts=(17,25,41))
    e=[x["relative_L2_error_curlA_axis"] for x in d["grid_convergence"]]
    assert e[-1] < e[0]
