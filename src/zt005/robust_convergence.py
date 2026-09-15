
from __future__ import annotations
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy

def sample_points(cyl,n,seed):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.82*cyl.radius_m,-0.82*cyl.height_m/2,-0.82*cyl.radius_m],
                      [0.82*cyl.radius_m,0.82*cyl.height_m/2,0.82*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.82*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def wire_distance(points,sources):
    d=np.full(len(points),np.inf)
    for src in sources:
        w=np.asarray(src[0],float)
        d=np.minimum(d,np.min(np.linalg.norm(points[:,None,:]-w[None,:,:],axis=2),axis=1))
    return d

def field_T(name,points,backpack,drive,h):
    sources=build_sources(name,backpack); omega=2*np.pi*drive.frequency_hz
    def A_fn(p): return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]
    A=np.vstack([A_fn(p) for p in points])
    B=np.vstack([curl_point_4(A_fn,p,h) for p in points])
    E=omega*A
    return E,B,em_stress_energy(E,B)["T"],sources

def robust_compare(X,Y):
    d=np.linalg.norm(X-Y,axis=1)
    ref=np.linalg.norm(Y,axis=1)
    rel=d/np.maximum(ref,1e-30)
    return {
        "mean":float(np.mean(rel)),
        "median":float(np.median(rel)),
        "p95":float(np.quantile(rel,0.95)),
        "max":float(np.max(rel))
    }
