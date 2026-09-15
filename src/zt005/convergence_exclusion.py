import numpy as np
from scipy.spatial import cKDTree

def sampled_polyline_points(polyline, samples_per_segment=2):
    q=np.asarray(polyline,float)
    if len(q)<2:
        return q
    pts=[]
    for i in range(len(q)):
        a=q[i]; b=q[(i+1)%len(q)]
        ts=np.linspace(0,1,samples_per_segment,endpoint=False)
        pts.append(a[None,:] + ts[:,None]*(b-a)[None,:])
    return np.vstack(pts)

def nearest_wire_distance(points, polylines, samples_per_segment=2, workers=1):
    p=np.asarray(points,float)
    clouds=[sampled_polyline_points(line,samples_per_segment) for line in polylines]
    cloud=np.vstack(clouds) if clouds else np.empty((0,3))
    if len(cloud)==0:
        return np.full(len(p),np.inf)
    tree=cKDTree(cloud)
    d,_=tree.query(p,k=1,workers=workers)
    return d

def exclusion_report(points, polylines, dx_m, factors=(2.,3.,4.), workers=1):
    d=nearest_wire_distance(points,polylines,samples_per_segment=2,workers=workers)
    rows=[]
    for k in factors:
        keep=d >= k*dx_m
        rows.append({
            "factor":float(k),
            "threshold_m":float(k*dx_m),
            "excluded_fraction":float((~keep).mean()),
            "kept_fraction":float(keep.mean()),
            "dmin_min_m":float(d.min()),
            "dmin_p01_m":float(np.quantile(d,.01)),
            "dmin_median_m":float(np.median(d))
        })
    return rows
