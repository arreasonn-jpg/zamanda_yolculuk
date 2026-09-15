import numpy as np
from .constants import MU0

def vector_potential(points, sources, peak_current_a):
    """A0 for I(t)=I0 cos(omega t), evaluated from filament segments."""
    p = np.asarray(points, dtype=float)
    A = np.zeros_like(p)
    for pts, seg in sources:
        mid = pts + 0.5*seg
        r = p[:, None, :] - mid[None, :, :]
        R = np.linalg.norm(r, axis=2)
        R = np.maximum(R, 1e-12)
        dlI = seg * peak_current_a
        A += (MU0/(4*np.pi)) * np.sum(dlI[None,:,:] / R[:,:,None], axis=1)
    return A

def biot_savart(points, sources, peak_current_a):
    p = np.asarray(points, dtype=float)
    B = np.zeros_like(p)
    for pts, seg in sources:
        mid = pts + 0.5*seg
        r = p[:, None, :] - mid[None, :, :]
        R3 = np.maximum(np.linalg.norm(r, axis=2)**3, 1e-36)
        dlI = seg * peak_current_a
        B += (MU0/(4*np.pi))*np.sum(
            np.cross(dlI[None,:,:], r, axis=2)/R3[:,:,None], axis=1
        )
    return B

def curl_on_regular_grid(V, axes):
    """Central-difference curl for V[nx,ny,nz,3]."""
    x,y,z = axes
    dx = float(x[1]-x[0]); dy=float(y[1]-y[0]); dz=float(z[1]-z[0])
    dVz_dy = np.gradient(V[...,2], dy, axis=1, edge_order=2)
    dVy_dz = np.gradient(V[...,1], dz, axis=2, edge_order=2)
    dVx_dz = np.gradient(V[...,0], dz, axis=2, edge_order=2)
    dVz_dx = np.gradient(V[...,2], dx, axis=0, edge_order=2)
    dVy_dx = np.gradient(V[...,1], dx, axis=0, edge_order=2)
    dVx_dy = np.gradient(V[...,0], dy, axis=1, edge_order=2)
    C = np.empty_like(V)
    C[...,0] = dVz_dy - dVy_dz
    C[...,1] = dVx_dz - dVz_dx
    C[...,2] = dVy_dx - dVx_dy
    return C

def divergence_on_regular_grid(V, axes):
    x,y,z=axes
    dx=float(x[1]-x[0]); dy=float(y[1]-y[0]); dz=float(z[1]-z[0])
    return (
        np.gradient(V[...,0],dx,axis=0,edge_order=2)
        +np.gradient(V[...,1],dy,axis=1,edge_order=2)
        +np.gradient(V[...,2],dz,axis=2,edge_order=2)
    )

def harmonic_fields(A0, axes, omega):
    """At t=pi/(2 omega): E=omega*A0, B=0 for the cosine-current baseline.
    We also return the amplitude B0=curl(A0) for Faraday validation."""
    B0 = curl_on_regular_grid(A0, axes)
    E90 = omega*A0
    return B0, E90
