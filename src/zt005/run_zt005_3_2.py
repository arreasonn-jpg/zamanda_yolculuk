
import json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .maxwell import vector_potential, biot_savart
from .pointwise_curl import curl_point_2, curl_point_4

OUT=Path("results/exploratory/zt005_3_2_results.json")

def sample_points(cyl, n_radial=6, n_z=7, n_phi=16, seed=0):
    rng=np.random.default_rng(seed)
    pts=[]
    for y in np.linspace(-0.9*cyl.height_m/2,0.9*cyl.height_m/2,n_z):
        for r in np.linspace(0,0.85*cyl.radius_m,n_radial):
            phis=[0.0] if r==0 else np.linspace(0,2*np.pi,n_phi,endpoint=False)
            for phi in phis:
                pts.append([r*np.cos(phi),y,r*np.sin(phi)])
    pts=np.asarray(pts,float)

    extra=[]
    while len(extra)<200:
        q=rng.uniform(
            [-0.85*cyl.radius_m,-0.9*cyl.height_m/2,-0.85*cyl.radius_m],
            [ 0.85*cyl.radius_m, 0.9*cyl.height_m/2, 0.85*cyl.radius_m]
        )
        if q[0]**2+q[2]**2 <= (0.85*cyl.radius_m)**2:
            extra.append(q)
    return np.vstack([pts,np.asarray(extra)])

def rel_l2(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))

def main():
    cyl=ActiveCylinder(); backpack=Backpack(); drive=Drive()
    points=sample_points(cyl)
    hs=[0.02,0.01,0.005]

    print("="*84)
    print("ZAMANDA YOLCULUK — ZT-005.3.2")
    print("SPARSE / POINTWISE FIELD CROSS-VALIDATION")
    print("="*84)
    print(f"Validation points: {len(points)}")
    print(f"h values: {hs}")

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,backpack)
        B0=biot_savart(points,sources,drive.peak_current_a)

        def A_fn(p):
            return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]

        rows=[]
        for h in hs:
            c2=[]; c4=[]; bs=[]
            for p,b in zip(points,B0):
                safe=True
                for src in sources:
                    wire=np.asarray(src[0])
                    if np.min(np.linalg.norm(wire-p[None,:],axis=1)) < 3*h:
                        safe=False; break
                if not safe:
                    continue
                c2.append(curl_point_2(A_fn,p,h))
                c4.append(curl_point_4(A_fn,p,h))
                bs.append(b)
            c2=np.asarray(c2); c4=np.asarray(c4); bs=np.asarray(bs)
            if len(bs)==0:
                rows.append({"h_m":h,"points":0})
                continue
            e2=rel_l2(c2,bs); e4=rel_l2(c4,bs)
            rows.append({
                "h_m":h,
                "points":int(len(bs)),
                "curl2_vs_BS_relL2":e2,
                "curl4_vs_BS_relL2":e4,
                "curl4_improvement_factor":e2/max(e4,1e-30)
            })
        results[name]={"validation_points_total":len(points),"results":rows}
        valid=[r for r in rows if r["points"]]
        for r in valid:
            print(f"{name} h={r['h_m']:.4f} points={r['points']} "
                  f"curl2={r['curl2_vs_BS_relL2']:.6e} "
                  f"curl4={r['curl4_vs_BS_relL2']:.6e} "
                  f"improvement={r['curl4_improvement_factor']:.3f}x")

    OUT.write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: sparse field cross-validation; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
