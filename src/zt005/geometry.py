
import numpy as np
from .model import Backpack

def _pack_point(b, x, y, z):
    return np.array([x, b.center_y_m + y, z], dtype=float)

def circle_loop(radius, center, axis="z", n=128):
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
    if axis == "z":
        pts = np.column_stack([radius*np.cos(t), radius*np.sin(t), np.zeros_like(t)])
    elif axis == "y":
        pts = np.column_stack([radius*np.cos(t), np.zeros_like(t), radius*np.sin(t)])
    elif axis == "x":
        pts = np.column_stack([np.zeros_like(t), radius*np.cos(t), radius*np.sin(t)])
    else:
        raise ValueError(axis)
    pts += np.asarray(center)
    seg = np.roll(pts, -1, axis=0) - pts
    return pts, seg

def helical_loop(radius, length, turns, center, axis="y", n=240):
    t = np.linspace(0, 2*np.pi*turns, n, endpoint=False)
    axial = (length/(2*np.pi*turns))*t - length/2
    if axis == "y":
        pts = np.column_stack([radius*np.cos(t), axial, radius*np.sin(t)])
    elif axis == "z":
        pts = np.column_stack([radius*np.cos(t), radius*np.sin(t), axial])
    else:
        raise ValueError(axis)
    pts += np.asarray(center)
    seg = np.roll(pts, -1, axis=0) - pts
    return pts, seg

def build_sources(name, b: Backpack, circle_n=128, helix_n=240):
    if name == "G1":
        return [circle_loop(0.065, _pack_point(b, 0, y, 0), "y", circle_n)
                for y in np.linspace(-0.055, 0.055, 6)]
    if name == "G2":
        return [circle_loop(0.065, _pack_point(b, x, 0, 0), "x", circle_n)
                for x in np.linspace(-0.12, 0.12, 4)]
    if name == "G3":
        a = helical_loop(0.065, 0.08, 1.0, _pack_point(b, -0.075, 0.0, 0), "y", helix_n)
        c = helical_loop(0.065, 0.08, 1.0, _pack_point(b,  0.075, 0.0, 0), "y", helix_n)
        return [(a[0], a[1]), (c[0], -c[1])]
    if name == "G4":
        return [
            circle_loop(0.065, _pack_point(b, 0, -0.055, -0.08), "y", circle_n),
            circle_loop(0.065, _pack_point(b, 0,  0.055,  0.08), "y", circle_n),
            circle_loop(0.065, _pack_point(b, 0, 0, 0), "x", circle_n),
        ]
    raise ValueError(name)

def check_source_envelope(sources, b: Backpack, tol=1e-10):
    lo = np.array([-b.width_m/2, b.center_y_m-b.depth_m/2, -b.height_m/2])
    hi = np.array([ b.width_m/2, b.center_y_m+b.depth_m/2,  b.height_m/2])
    for pts, _ in sources:
        if np.any(pts < lo-tol) or np.any(pts > hi+tol):
            return False
    return True
