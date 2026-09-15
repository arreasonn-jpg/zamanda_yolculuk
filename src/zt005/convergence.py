
import numpy as np
from .constants import MU0
from .geometry import circle_loop
from .maxwell import vector_potential, biot_savart

def analytic_loop_axis_B(z, radius, current):
    z = np.asarray(z, dtype=float)
    # circle_loop(axis="y") is parameterized with +t orientation; its
    # magnetic field points in -y at the origin by the right-hand rule.
    return -MU0 * current * radius**2 / (2.0 * (radius**2 + z**2)**1.5)

def axis_curlA_y(zvals, radius, current, nseg, h):
    """curl(A)_y = dA_x/dz - dA_z/dx on the symmetry axis.

    Uses a local central-difference stencil at x=z_cross=0.  Only a tiny
    number of A evaluations are required, avoiding an unnecessarily huge
    3-D mesh for the convergence benchmark.
    """
    zvals=np.asarray(zvals,dtype=float)
    sources=[circle_loop(radius, np.zeros(3), "y", nseg)]
    p_xp=np.column_stack([np.full_like(zvals,h),zvals,np.zeros_like(zvals)])
    p_xm=np.column_stack([np.full_like(zvals,-h),zvals,np.zeros_like(zvals)])
    p_zp=np.column_stack([np.zeros_like(zvals),zvals,np.full_like(zvals,h)])
    p_zm=np.column_stack([np.zeros_like(zvals),zvals,np.full_like(zvals,-h)])
    A_xp=vector_potential(p_xp,sources,current)
    A_xm=vector_potential(p_xm,sources,current)
    A_zp=vector_potential(p_zp,sources,current)
    A_zm=vector_potential(p_zm,sources,current)
    dAx_dz=(A_zp[:,0]-A_zm[:,0])/(2*h)
    dAz_dx=(A_xp[:,2]-A_xm[:,2])/(2*h)
    return dAx_dz-dAz_dx

def benchmark_single_loop(segment_counts=(32,64,128,256,512),
                           grid_counts=(17,25,41,61),
                           radius=0.12, current=100.0):
    out={"segment_convergence":[], "grid_convergence":[]}

    zvals=np.array([-0.30,-0.18,-0.06,0.06,0.18,0.30])
    points=np.column_stack([np.zeros_like(zvals),zvals,np.zeros_like(zvals)])
    exact=analytic_loop_axis_B(zvals,radius,current)

    for nseg in segment_counts:
        sources=[circle_loop(radius, np.zeros(3), "y", nseg)]
        B=biot_savart(points,sources,current)
        rel=np.linalg.norm(B[:,1]-exact)/np.linalg.norm(exact)
        out["segment_convergence"].append({
            "segments":int(nseg),
            "relative_L2_error":float(rel)
        })

    # Map grid resolution to the local finite-difference spacing.
    # The physical benchmark interval is fixed at [-0.18,0.18].
    for n in grid_counts:
        h=0.36/(n-1)
        pred=axis_curlA_y(zvals,radius,current,512,h)
        rel=np.linalg.norm(pred-exact)/np.linalg.norm(exact)
        out["grid_convergence"].append({
            "grid_n":int(n),
            "grid_spacing_m":float(h),
            "relative_L2_error_curlA_axis":float(rel)
        })
    return out
