
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .batch_em_source import compute_T_batch
from .cyl_quadrature import cylindrical_midpoints,quadrature_integral_many,EINSTEIN_GREEN_FACTOR
from .spatial_gr import inverse_trace_reverse

OUT=Path("results/exploratory/zt006_2_2_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)
    reslist=[(4,8,4),(6,12,6),(8,16,8)] if a.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16)]
    print("="*92)
    print("ZAMANDA YOLCULUK — ZT-006.2.2")
    print("PHYSICS-CORRECTED + BATCHED DETERMINISTIC QUADRATURE")
    print("="*92)
    print(f"Mode: {a.mode}")
    print(f"Einstein Green factor 4G/c^4 = {EINSTEIN_GREEN_FACTOR:.6e}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}]")
        rows=[]
        for nr,np_,nz in reslist:
            pts,w=cylindrical_midpoints(cyl.radius_m,cyl.height_m,nr,np_,nz)
            T=compute_T_batch(name,pts,b,d,h=.005)
            hb=quadrature_integral_many(obs,pts,T,w,softening_m=.0025)
            hp=np.stack([inverse_trace_reverse(x) for x in hb])
            h00=hp[:,0,0]; h0=np.linalg.norm(hp[:,0,1:],axis=1)
            rows.append({"resolution":[nr,np_,nz],"source_points":len(pts),
                         "center_h00":float(h00[0]),"center_h0i_norm":float(h0[0]),
                         "h00":h00.tolist(),"h0i_norm":h0.tolist()})
            print(f"  ({nr},{np_},{nz}) N={len(pts)} center h00={h00[0]:.6e} |h0i|={h0[0]:.6e}")
        for i in range(1,len(rows)):
            p=rows[i-1]["center_h00"]; c=rows[i]["center_h00"]
            rows[i]["relative_change_center_h00"]=abs(c-p)/max(abs(c),1e-300)
        rows[0]["relative_change_center_h00"]=None
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.2","mode":a.mode,"observers":obs.tolist(),"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: corrected deterministic source integration; no Earth/CTC conclusion.")

if __name__=="__main__": main()
