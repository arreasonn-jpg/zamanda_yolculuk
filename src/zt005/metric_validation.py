import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])

def metric_from_hbar(hbar):
    h = np.asarray(hbar, dtype=float)
    if h.shape != (4, 4):
        raise ValueError("hbar must be 4x4")
    return ETA + 0.5 * (h + h.T)

def metric_invariants(g):
    g = 0.5 * (np.asarray(g, dtype=float) + np.asarray(g, dtype=float).T)
    eig = np.linalg.eigvalsh(g)
    return {
        "g00": float(g[0, 0]),
        "spatial_trace": float(np.trace(g[1:, 1:])),
        "min_eigenvalue": float(eig.min()),
        "max_eigenvalue": float(eig.max()),
        "det": float(np.linalg.det(g)),
        "negative_eigenvalues": int(np.sum(eig < 0.0)),
        "positive_eigenvalues": int(np.sum(eig > 0.0)),
        "lorentzian_signature": bool(np.sum(eig < 0.0) == 1 and np.sum(eig > 0.0) == 3),
        "det_distance_from_minkowski": float(abs(np.linalg.det(g) + 1.0)),
    }

def timelike_check(g, v=None):
    if v is None:
        v = np.array([1.0, 0.0, 0.0, 0.0])
    else:
        v = np.asarray(v, dtype=float)
    interval = float(v @ g @ v)
    return {"interval": interval, "timelike": interval < 0.0}

def diagnostics(hbar):
    h = np.asarray(hbar, dtype=float)
    g = metric_from_hbar(h)
    return {
        "h_frobenius": float(np.linalg.norm(h)),
        "h_max_abs": float(np.max(np.abs(h))),
        "metric_invariants": metric_invariants(g),
        "rest_frame_timelike": timelike_check(g),
    }
