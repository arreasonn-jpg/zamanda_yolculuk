import numpy as np
from zt005.model import ActiveCylinder, Backpack, Drive
from zt005.gauss_cyl_quadrature import gauss_cylindrical_nodes
from zt005.green_convergence import focused_source_nodes, fast_regularized_em_fields
from zt005.green_solvers import solver_a_tensor_product
from zt005.geometry import build_sources

cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
srcs = build_sources("G1", bk)
print("cyl R,h =", cyl.radius_m, cyl.height_m)
ap = np.vstack([s[0] for s in srcs])
print("src x", ap[:,0].min(), ap[:,0].max())
print("src y", ap[:,1].min(), ap[:,1].max())
print("src z", ap[:,2].min(), ap[:,2].max())

obs = np.array([[0.,0.,2.],[0.,0.,3.],[0.,0.,5.],[2.,0.,0.],[3.,0.,0.]])

for label, fn in [("focused", focused_source_nodes), ("full", gauss_cylindrical_nodes)]:
    print("\n--", label, "--")
    prev = None
    for n in [6, 8, 12, 16, 20, 24]:
        pts, w = fn(cyl.radius_m, cyl.height_m, n, 2*n, n)
        T = fast_regularized_em_fields(pts, srcs, dv.peak_current_a, dv.frequency_hz)
        h = solver_a_tensor_product(obs, pts, w, T)
        e = None if prev is None else float(np.linalg.norm(h-prev)/max(np.linalg.norm(h),1e-300))
        print(f"n={n:3d} N={len(pts):6d} |h|={np.linalg.norm(h):.3e} rel={e}")
        prev = h
