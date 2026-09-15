
from __future__ import annotations
import numpy as np

def gauss_disk_nodes(nr, nphi):
    x,w=np.polynomial.legendre.leggauss(nr)
    r=(x+1)/2
    wr=w/2
    phi=2*np.pi*(np.arange(nphi)+0.5)/nphi
    wphi=2*np.pi/nphi
    pts=[]; weights=[]
    for i,ri in enumerate(r):
        for ph in phi:
            pts.append((ri*np.cos(ph),ri*np.sin(ph)))
            weights.append(2*ri*wr[i]*(wphi/(2*np.pi)))
    weights=np.asarray(weights,float)
    return np.asarray(pts,float),weights/weights.sum()

def _tube_basis(tangent, preferred_axis):
    t=np.asarray(tangent,float); n=np.linalg.norm(t)
    if n<1e-15: t=np.array([1.,0.,0.]); n=1.
    t=t/n
    a=np.asarray(preferred_axis,float)
    if abs(np.dot(t,a))>0.9:
        a=np.array([0.,0.,1.])
        if abs(np.dot(t,a))>0.9: a=np.array([1.,0.,0.])
    nv=np.cross(t,a); nv/=max(np.linalg.norm(nv),1e-15)
    bv=np.cross(t,nv); bv/=max(np.linalg.norm(bv),1e-15)
    return nv,bv

def build_area_weighted_groups(name,backpack,radius_m,nr=2,nphi=8,circle_n=128,helix_n=240):
    from .geometry import build_sources
    base=build_sources(name,backpack,circle_n=circle_n,helix_n=helix_n)
    groups=[]
    if radius_m<=0:
        return [[(src,1.0)] for src in base]
    offsets,weights=gauss_disk_nodes(nr,nphi)
    for src in base:
        pts,seg=src; pts=np.asarray(pts,float); seg=np.asarray(seg,float)
        group=[]
        for (aa,bb),wt in zip(offsets,weights):
            new=np.empty_like(pts)
            for k in range(len(pts)):
                nv,bv=_tube_basis(seg[k-1]+seg[k],np.array([0.,1.,0.]))
                new[k]=pts[k]+radius_m*(aa*nv+bb*bv)
            group.append(((new,seg.copy()),float(wt)))
        groups.append(group)
    return groups
