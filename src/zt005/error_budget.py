
from __future__ import annotations
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy
from .cyl_quadrature import cylindrical_midpoints, quadrature_integral_many
from .spatial_gr import inverse_trace_reverse

def validation_points(radius=0.75, height=2.0, nr=4, nphi=12, nz=5):
    pts=[]
    for y in np.linspace(-height/2,height/2,nz):
        for r in np.linspace(0,radius,nr):
            phis=[0.] if r==0 else np.linspace(0,2*np.pi,nphi,endpoint=False)
            for p in phis:
                pts.append([r*np.cos(p),y,r*np.sin(p)])
    return np.asarray(pts,float)

def T_at_fixed_points(name, points, backpack, drive, h):
    sources=build_sources(name,backpack)
    omega=2*np.pi*drive.frequency_hz
    def A_fn(p):
        return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]
    A=np.vstack([A_fn(p) for p in points])
    B=np.vstack([curl_point_4(A_fn,p,h) for p in points])
    E=omega*A
    return em_stress_energy(E,B)["T"]

def relative_tensor_l2(a,b):
    return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))

def softening_sensitivity(obs, source_points, T, weights, softenings):
    out=[]
    for s in softenings:
        hb=quadrature_integral_many(obs,source_points,T,weights,softening_m=s)
        hp=np.stack([inverse_trace_reverse(x) for x in hb])
        out.append({
            "softening_m":float(s),
            "center_h00":float(hp[0,0,0]),
            "center_h0i_norm":float(np.linalg.norm(hp[0,0,1:]))
        })
    ref=out[-1]
    for row in out:
        row["relative_change_from_last_h00"]=abs(row["center_h00"]-ref["center_h00"])/max(abs(ref["center_h00"]),1e-300)
    return out

def kernel_resolution_with_frozen_T(name, cyl, backpack, drive, h_field, resolutions, obs):
    # For each resolution, evaluate T on its own quadrature grid using the SAME
    # field-reconstruction h. The diagnostic compares raw Green contribution and
    # full h. This is a practical separation report, not a mathematically exact
    # operator decomposition.
    rows=[]
    for nr,nphi,nz in resolutions:
        pts,w=cylindrical_midpoints(cyl.radius_m,cyl.height_m,nr,nphi,nz)
        T=T_at_fixed_points(name,pts,backpack,drive,h_field)
        hb=quadrature_integral_many(obs,pts,T,w,softening_m=0.0025)
        hp=np.stack([inverse_trace_reverse(x) for x in hb])
        rows.append({
            "resolution":[nr,nphi,nz],
            "points":int(len(pts)),
            "center_h00":float(hp[0,0,0]),
            "center_h0i_norm":float(np.linalg.norm(hp[0,0,1:])),
        })
    for i in range(1,len(rows)):
        a=rows[i-1]["center_h00"]; b=rows[i]["center_h00"]
        rows[i]["relative_change_h00"]=abs(b-a)/max(abs(b),1e-300)
    rows[0]["relative_change_h00"]=None
    return rows
