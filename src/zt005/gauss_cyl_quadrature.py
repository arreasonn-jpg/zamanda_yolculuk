
from __future__ import annotations
import numpy as np

C = 299_792_458.0
G = 6.67430e-11
EINSTEIN_GREEN_FACTOR = 4.0 * G / C**4

def gauss_cylindrical_nodes(radius_m, height_m, nr, nphi, nz):
    # Gauss-Legendre in r and z; uniform trapezoidal in periodic phi.
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xz, wz = np.polynomial.legendre.leggauss(nz)

    r = 0.5*radius_m*(xr + 1.0)
    wr = 0.5*radius_m*wr

    z = 0.5*height_m*xz
    wz = 0.5*height_m*wz

    phi = 2*np.pi*(np.arange(nphi)+0.5)/nphi
    wphi = 2*np.pi/nphi

    rr, pp, zz = np.meshgrid(r,phi,z,indexing="ij")
    pts=np.column_stack([
        (rr*np.cos(pp)).ravel(),
        zz.ravel(),
        (rr*np.sin(pp)).ravel()
    ])

    weights=(wr[:,None,None]*wphi*wz[None,None,:]*rr).ravel()
    return pts, weights

def integrate_green(obs, source_points, source_values, weights, softening_m=0.0):
    obs=np.asarray(obs,float)
    src=np.asarray(source_points,float)
    vals=np.asarray(source_values,float)
    w=np.asarray(weights,float)
    d=np.linalg.norm(src-obs[None,:],axis=1)
    if softening_m>0:
        d=np.sqrt(d*d+softening_m**2)
    d=np.maximum(d,1e-15)
    raw = np.sum((w/d)[:,None,None]*vals,axis=0)
    return EINSTEIN_GREEN_FACTOR * raw

def integrate_green_many(obs_points, source_points, source_values, weights, softening_m=0.0):
    return np.stack([integrate_green(p,source_points,source_values,weights,softening_m)
                     for p in np.asarray(obs_points)])
