
import json,argparse
from pathlib import Path
import numpy as np

from .model import ActiveCylinder,Backpack,Drive
from .batch_em_source import compute_T_batch
from .cyl_quadrature import cylindrical_midpoints
from .gauss_cyl_quadrature import gauss_cylindrical_nodes,integrate_green_many,EINSTEIN_GREEN_FACTOR
from .spatial_gr import inverse_trace_reverse

OUT=Path("results/exploratory/zt006_2_9_1_results.json")

def metric(obs,pts,w,T):
    hbar=integrate_green_many(obs,pts,T,w,softening_m=0.0025)
    return np.stack([inverse_trace_reverse(x) for x in hbar])

def as_serializable_metric(m):
    m=np.asarray(m)
    return {
        "center_h00":float(m[0,0,0]),
        "center_h0i_norm":float(np.linalg.norm(m[0,0,1:])),
        "all_h00":[float(x) for x in m[:,0,0]],
        "all_h0i_norm":[float(np.linalg.norm(x[0,1:])) for x in m],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)
    res=[(4,8,4),(6,12,6),(8,16,8)] if a.mode=="quick" else [(6,12,6),(8,16,8),(12,24,12),(16,32,16)]
    print("="*100)
    print(f"Einstein Green factor 4G/c^4 = {EINSTEIN_GREEN_FACTOR:.6e}")
    print("ZAMANDA YOLCULUK — ZT-006.2.9.1")
    print("GAUSS-LEGENDRE / PERIODIC PHI GREEN-QUADRATURE + JSON-SAFE OUTPUT")
    print("="*100)

    results={}
    for name in ["G1","G2","G3","G4"]:
        rows=[]
        print(f"\n[{name}]")
        for nr,nphi,nz in res:
            pts,w=gauss_cylindrical_nodes(cyl.radius_m,cyl.height_m,nr,nphi,nz)
            T=compute_T_batch(name,pts,b,d,h=.005)
            m=metric(obs,pts,w,T)
            data=as_serializable_metric(m)
            data["resolution"]=[nr,nphi,nz]
            data["source_points"]=int(len(pts))
            rows.append(data)
            print(f"  Gauss ({nr},{nphi},{nz}) N={len(pts)} "
                  f"center h00={data['center_h00']:.6e} "
                  f"|h0i|={data['center_h0i_norm']:.6e}")

        ref=np.array(rows[-1]["all_h00"],float)
        for row in rows:
            cur=np.array(row["all_h00"],float)
            rel=np.linalg.norm(cur-ref)/max(np.linalg.norm(ref),1e-300)
            row["relative_h00_vector_error_vs_finest"]=float(rel)

        # Midpoint comparison at matched-ish node count is useful as a bias diagnostic.
        mid_res=res[min(2,len(res)-1)]
        ptsm,wm=cylindrical_midpoints(cyl.radius_m,cyl.height_m,*mid_res)
        Tm=compute_T_batch(name,ptsm,b,d,h=.005)
        mm=metric(obs,ptsm,wm,Tm)
        rows.append({"method":"midpoint_reference",
                     "resolution":list(mid_res),
                     **as_serializable_metric(mm)})
        results[name]=rows

    OUT.write_text(json.dumps({
        "stage":"ZT-006.2.9.1",
        "mode":a.mode,
        "observers":obs.tolist(),
        "results":results
    },indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: deterministic quadrature-method comparison; no GR/CTC conclusion.")

if __name__=="__main__": main()
