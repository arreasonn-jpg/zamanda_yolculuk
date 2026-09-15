
import numpy as np
from .constants import MU0
from .maxwell import vector_potential, biot_savart, curl_on_regular_grid
from .geometry import build_sources, check_source_envelope
from .metrics import rel_l2, rms
from .model import Backpack

def grid(cyl,nx,ny,nz):
    x=np.linspace(-cyl.radius_m,cyl.radius_m,nx)
    y=np.linspace(-cyl.height_m/2,cyl.height_m/2,ny)
    z=np.linspace(-cyl.radius_m,cyl.radius_m,nz)
    X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    mask=X*X+Z*Z <= cyl.radius_m**2
    pts=np.column_stack([X.ravel(),Y.ravel(),Z.ravel()])
    return x,y,z,X,Y,Z,pts,mask

def _chunked_field(func, points, sources, current, chunk=4096):
    out=np.empty((len(points),3),dtype=float)
    for i in range(0,len(points),chunk):
        out[i:i+chunk]=func(points[i:i+chunk],sources,current)
    return out

def evaluate_geometry(name,cyl,b,current,grid_shape=(25,17,25),
                      circle_n=128,helix_n=240,chunk=4096):
    nx,ny,nz=grid_shape
    x,y,z,X,Y,Z,pts,mask=grid(cyl,nx,ny,nz)
    sources=build_sources(name,b,circle_n,helix_n)
    if not check_source_envelope(sources,b):
        raise RuntimeError(f"{name} violates hard backpack envelope")
    Aflat=_chunked_field(vector_potential,pts,sources,current,chunk)
    Bbsflat=_chunked_field(biot_savart,pts,sources,current,chunk)
    A=Aflat.reshape(X.shape+(3,))
    Bbs=Bbsflat.reshape(X.shape+(3,))
    Bcurl=curl_on_regular_grid(A,(x,y,z))
    interior=np.zeros(X.shape,dtype=bool)
    interior[1:-1,1:-1,1:-1]=True
    interior &= mask
    # field magnitudes
    Bm=np.linalg.norm(Bcurl,axis=-1)
    Bbs_m=np.linalg.norm(Bbs,axis=-1)
    # normalized divergence metric
    dx=min(abs(x[1]-x[0]),abs(y[1]-y[0]),abs(z[1]-z[0]))
    dBx=np.gradient(Bcurl[...,0],x,axis=0,edge_order=2)
    dBy=np.gradient(Bcurl[...,1],y,axis=1,edge_order=2)
    dBz=np.gradient(Bcurl[...,2],z,axis=2,edge_order=2)
    div=dBx+dBy+dBz
    div_rel=rms(div,interior)/(max(rms(Bcurl,interior)/dx,1e-30))
    return {
        "geometry":name,
        "grid":list(map(int,grid_shape)),
        "circle_segments":int(circle_n),
        "helix_segments":int(helix_n),
        "points":int(len(pts)),
        "B_biot_savart_mean_T":float(Bbs_m[interior].mean()),
        "B_curlA_mean_T":float(Bm[interior].mean()),
        "B_curlA_vs_BS_relL2":rel_l2(Bcurl,Bbs,interior),
        "divB_relative":float(div_rel),
    }
