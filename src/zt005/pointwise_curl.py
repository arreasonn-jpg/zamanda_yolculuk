
import numpy as np

def axis_stencil_points(center, h, axis):
    c = np.asarray(center, dtype=float)
    p2 = c.copy(); p2[axis] += 2*h
    p1 = c.copy(); p1[axis] += h
    m1 = c.copy(); m1[axis] -= h
    m2 = c.copy(); m2[axis] -= 2*h
    return p2, p1, m1, m2

def curl_point_2(field_fn, point, h):
    p = np.asarray(point, dtype=float)
    deriv = []
    for axis in range(3):
        _,p1,m1,_ = axis_stencil_points(p,h,axis)
        deriv.append((field_fn(p1)-field_fn(m1))/(2*h))
    d_dx,d_dy,d_dz = deriv
    return np.array([
        d_dy[2] - d_dz[1],
        d_dz[0] - d_dx[2],
        d_dx[1] - d_dy[0],
    ], dtype=float)

def curl_point_4(field_fn, point, h):
    p = np.asarray(point, dtype=float)
    deriv = []
    for axis in range(3):
        p2,p1,m1,m2 = axis_stencil_points(p,h,axis)
        deriv.append(
            (-field_fn(p2) + 8*field_fn(p1)
             - 8*field_fn(m1) + field_fn(m2))/(12*h)
        )
    d_dx,d_dy,d_dz = deriv
    return np.array([
        d_dy[2] - d_dz[1],
        d_dz[0] - d_dx[2],
        d_dx[1] - d_dy[0],
    ], dtype=float)
