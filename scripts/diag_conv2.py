import numpy as np
from zt005.model import ActiveCylinder, Backpack, Drive
from zt005.gauss_cyl_quadrature import gauss_cylindrical_nodes
from zt005.green_convergence import focused_source_nodes, fast_regularized_em_fields
from zt005.green_solvers import solver_a_tensor_product
from zt005.geometry import build_sources
from zt005.manufactured_sources import GaussianBump, conserved_T

cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
obs = np.array([[0.,0.,2.],[0.,0.,3.],[0.,0.,5.],[2.,0.,0.],[3.,0.,0.]])

print("== TEST A: smooth Gaussian source (solver check) ==")
bump = GaussianBump(A=1.0, sigma_t=1e-3, sigma_x=0.5)
prev=None
for n in [8,12,16,20,24,32]:
    pts,w = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, n, 2*n, n)
    T = conserved_T(bump, np.zeros(len(pts)), pts)
    h = solver_a_tensor_product(obs, pts, w, T)
    e = None if prev is None else float(np.linalg.norm(h-prev)/max(np.linalg.norm(h),1e-300))
    print(f"n={n:3d} N={len(pts):6d} |h|={np.linalg.norm(h):.3e} rel={e}")
    prev=h

print("\n== TEST B: device filament, a_wire sweep ==")
srcs = build_sources("G1", bk)
for aw in [0.001, 0.005, 0.02, 0.05, 0.1]:
    print(f" a_wire={aw}")
    prev=None
    for n in [8,12,16,20]:
        pts,w = focused_source_nodes(cyl.radius_m, cyl.height_m, n, 2*n, n)
        T = fast_regularized_em_fields(pts, srcs, dv.peak_current_a, dv.frequency_hz, a_wire=aw)
        h = solver_a_tensor_product(obs, pts, w, T)
        e = None if prev is None else float(np.linalg.norm(h-prev)/max(np.linalg.norm(h),1e-300))
        print(f"  n={n:3d} |h|={np.linalg.norm(h):.3e} rel={e}")
        prev=h
