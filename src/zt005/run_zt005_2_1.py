import json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder, Backpack
from .geometry import build_sources
from .convergence_exclusion import exclusion_report

def make_grid(cyl,n):
    x=np.linspace(-cyl.radius_m,cyl.radius_m,n)
    y=np.linspace(-cyl.height_m/2,cyl.height_m/2,n+10)
    z=np.linspace(-cyl.radius_m,cyl.radius_m,n)
    X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    mask=(X*X+Z*Z)<=cyl.radius_m**2
    pts=np.column_stack((X[mask],Y[mask],Z[mask]))
    dx=min(x[1]-x[0],y[1]-y[0],z[1]-z[0])
    return pts,dx

def main():
    cyl=ActiveCylinder(); b=Backpack()
    results={"stage":"ZT-005.2.1-optimized","grids":[25,41,61],"geometries":{}}
    print("="*80)
    print("ZAMANDA YOLCULUK — ZT-005.2.1")
    print("OPTIMIZED WIRE-DISTANCE / INDEPENDENT-GRID DIAGNOSTICS")
    print("="*80)
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        lines=[np.asarray(s[0]) for s in sources]
        rows=[]
        for n in [25,41,61]:
            pts,dx=make_grid(cyl,n)
            rep=exclusion_report(pts,lines,dx,workers=1)
            rows.append({"grid":n,"points":int(len(pts)),"dx_m":float(dx),"exclusion":rep})
            k3=rep[1]
            print(f"{name} grid={n:>2} points={len(pts):>6} "
                  f"dx={dx:.4e} m excluded(k=3)={100*k3['excluded_fraction']:.2f}% "
                  f"p01={k3['dmin_p01_m']:.4e} m")
        results["geometries"][name]=rows
    Path("zt005_2_1_results.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(f"Results written to zt005_2_1_results.json")
    print("SCIENTIFIC STATUS: diagnostic only; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
