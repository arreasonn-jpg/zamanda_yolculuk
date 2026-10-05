
import json,argparse
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .geometry import build_sources
from .batch_em_source import compute_T_batch
from .cyl_quadrature import cylindrical_midpoints
from .frozen_source_gr import frozen_green_metric,relative_metric_error

OUT=Path("results/exploratory/zt006_2_9_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0.,0.,0.],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)
    grids=[(6,12,6),(8,16,8),(12,24,12)] if a.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16),(20,40,20)]
    print("="*100)
    print("ZAMANDA YOLCULUK — ZT-006.2.9")
    print("FROZEN T_MUNU + DETERMINISTIC GREEN-INTEGRAL CONVERGENCE")
    print("="*100)
    print(f"Grids: {grids}")
    print("Purpose: isolate Green quadrature convergence from random source sampling.")
    results={}
    for name in ["G1","G2","G3","G4"]:
        rows=[]
        for res in grids:
            nr,nphi,nz=res
            pts,w=cylindrical_midpoints(cyl.radius_m,cyl.height_m,nr,nphi,nz)
            T=compute_T_batch(name,pts,b,d,h=0.005)
            metric=frozen_green_metric(obs,pts,w,T,softening_m=0.0025)
            rows.append({"resolution":list(res),"source_points":len(pts),"metric":metric})
            print(f"{name} grid={res} N={len(pts)} center h00={metric[0,0,0]:.6e} |h0i|={np.linalg.norm(metric[0,0,1:]):.6e}")
        ref=rows[-1]["metric"]
        for row in rows:
            row["error_vs_finest"]=relative_metric_error(row["metric"],ref)
        results[name]=rows
        print(f"{name} median errors vs finest: {[r['error_vs_finest']['median'] for r in rows]}")
    OUT.write_text(json.dumps({"stage":"ZT-006.2.9","mode":a.mode,"observers":obs.tolist(),"results":results},indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: frozen-source-grid separation study; no GR/CTC conclusion.")

if __name__=="__main__": main()
