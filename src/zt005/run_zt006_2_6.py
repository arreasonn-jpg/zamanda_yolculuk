
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .local_conditioning import field_values,local_conditioning

OUT=Path("zt006_2_6_results.json")

def sample_points(cyl,n=300,seed=321):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,120 if args.mode=="quick" else 300)
    hs=[0.02,0.01,0.005,0.0025]
    print("="*96)
    print("ZAMANDA YOLCULUK — ZT-006.2.6")
    print("FIELD GRADIENT / LOCAL CONDITIONING STUDY")
    print("="*96)
    print(f"Validation points: {len(pts)}")
    print(f"h values: {hs}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        rows=[]
        for h in hs:
            E,B=field_values(name,pts,b,d,h)
            c=local_conditioning(E,B,pts)
            rows.append({"h_m":h,**c})
            print(f"{name} h={h:.4f} gradB/B={c['gradB_over_B_mean']:.6e} "
                  f"p95={c['gradB_over_B_p95']:.6e} gradE/E={c['gradE_over_E_mean']:.6e}")
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.6","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: local field-conditioning diagnostics; no GR/CTC conclusion.")

if __name__=="__main__": main()
