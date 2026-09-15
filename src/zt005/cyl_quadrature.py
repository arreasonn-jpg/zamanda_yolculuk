
from __future__ import annotations
import numpy as np

C = 299_792_458.0
G = 6.67430e-11
EINSTEIN_GREEN_FACTOR = 4.0 * G / C**4

def cylindrical_midpoints(radius_m, height_m, nr, nphi, nz):
    # Midpoint rule in r, phi, z. The radial Jacobian r is returned separately.
    dr = radius_m / nr
    dphi = 2.0*np.pi / nphi
    dz = height_m / nz

    r = (np.arange(nr) + 0.5) * dr
    phi = (np.arange(nphi) + 0.5) * dphi
    z = -height_m/2 + (np.arange(nz) + 0.5) * dz

    rr, pp, zz = np.meshgrid(r, phi, z, indexing="ij")
    x = rr*np.cos(pp)
    y = zz
    zc = rr*np.sin(pp)

    pts = np.column_stack((x.ravel(), y.ravel(), zc.ravel()))
    weights = (rr * dr * dphi * dz).ravel()
    return pts, weights

def quadrature_integral(obs, source_pts, source_values, weights,
                        softening_m=0.0):
    obs = np.asarray(obs, dtype=float)
    src = np.asarray(source_pts, dtype=float)
    vals = np.asarray(source_values, dtype=float)
    w = np.asarray(weights, dtype=float)

    d = np.linalg.norm(src - obs[None, :], axis=1)
    if softening_m > 0:
        d = np.sqrt(d*d + softening_m*softening_m)
    d = np.maximum(d, 1e-15)

    integral = np.sum((w/d)[:, None, None] * vals, axis=0)
    return EINSTEIN_GREEN_FACTOR * integral

def quadrature_integral_many(obs_points, source_pts, source_values, weights,
                             softening_m=0.0):
    return np.stack([
        quadrature_integral(p, source_pts, source_values, weights, softening_m)
        for p in np.asarray(obs_points)
    ])
