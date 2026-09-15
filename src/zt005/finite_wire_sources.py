
from __future__ import annotations
import numpy as np

def _tube_basis(tangent, preferred_axis):
    t=np.asarray(tangent,float); nrm=np.linalg.norm(t)
    if nrm<1e-15: t=np.array([1.,0.,0.]); nrm=1.
    t=t/nrm; a=np.asarray(preferred_axis,float)
    if abs(np.dot(t,a))>0.9:
        a=np.array([0.,0.,1.])
        if abs(np.dot(t,a))>0.9: a=np.array([1.,0.,0.])
    n=np.cross(t,a); n/=max(np.linalg.norm(n),1e-15)
    b=np.cross(t,n); b/=max(np.linalg.norm(b),1e-15)
    return n,b

def _offset_samples(radius_m,samples):
    if radius_m<=0 or samples<=1: return [(0.,0.)]
    out=[(0.,0.)]; remaining=samples-1; ring=1
    rings=max(1,int(np.ceil((samples-1)/6)))
    while remaining>0:
        m=min(6*ring,remaining)
        rr=radius_m*ring/rings*0.85
        for j in range(m):
            ang=2*np.pi*j/m
            out.append((rr*np.cos(ang),rr*np.sin(ang)))
        remaining-=m; ring+=1
    return out

def thicken_source(source,radius_m,samples=6):
    pts,seg=source; pts=np.asarray(pts,float); seg=np.asarray(seg,float)
    out=[]
    for aa,bb in _offset_samples(radius_m,samples):
        new=np.empty_like(pts)
        for k in range(len(pts)):
            n,b=_tube_basis(seg[k-1]+seg[k],np.array([0.,1.,0.]))
            new[k]=pts[k]+aa*n+bb*b
        out.append((new,seg.copy()))
    return out

def build_finite_sources_grouped(name,backpack,radius_m,cross_samples=6,circle_n=128,helix_n=240):
    from .geometry import build_sources
    return [thicken_source(src,radius_m,cross_samples)
            for src in build_sources(name,backpack,circle_n=circle_n,helix_n=helix_n)]

def build_finite_sources(name,backpack,radius_m,cross_samples=6,circle_n=128,helix_n=240):
    return [s for g in build_finite_sources_grouped(name,backpack,radius_m,cross_samples,circle_n,helix_n) for s in g]
