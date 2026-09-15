
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .area_weighted_finite_wire import build_area_weighted_groups

OUT=Path("zt006_2_15_results.json")

def sample_points(cyl,n=100,seed=515):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def rel(a,b):
    return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))

def fields(name,points,b,d,radius,nr,nphi,h=.005):
    from .maxwell import vector_potential
    from .pointwise_curl import curl_point_4
    from .validated_stress_energy import em_stress_energy
    groups=build_area_weighted_groups(name,b,radius,nr,nphi)
    sources=[s for g in groups for s,w in g]
    currents=[d.peak_current_a*w for g in groups for s,w in g]
    def A_fn(p):
        p=np.atleast_2d(p); A=np.zeros((len(p),3))
        for src,I in zip(sources,currents): A+=vector_potential(p,[src],I)
        return A[0]
    p=np.asarray(points,float)
    A=np.vstack([A_fn(x) for x in p])
    B=np.vstack([curl_point_4(A_fn,x,h) for x in p])
    E=2*np.pi*d.frequency_hz*A
    T=em_stress_energy(E,B)["T"]
    return E,B,T,len(groups),len(sources)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,100 if a.mode=="quick" else 200)
    radii=[.001,.0025,.005] if a.mode=="quick" else [.001,.0025,.005,.0075]
    rules=[(1,6),(2,6),(3,8)] if a.mode=="quick" else [(1,6),(2,8),(3,10),(4,12)]
    results={}
    print("="*108)
    print("ZAMANDA YOLCULUK — ZT-006.2.15")
    print("AREA-WEIGHTED FINITE-CROSS-SECTION CONVERGENCE")
    print("="*108)
    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}]"); results[name]={}
        for r in radii:
            ref=fields(name,pts,b,d,r,*rules[-1])
            Eref,Bref,Tref=ref[:3]
            rows=[]
            for nr,nphi in rules:
                E,B,T,pc,sc=fields(name,pts,b,d,r,nr,nphi)
                row={"nr":nr,"nphi":nphi,"parent_sources":pc,"subfilaments":sc,
                     "E_rel_vs_ref":rel(E,Eref),"B_rel_vs_ref":rel(B,Bref),
                     "T_rel_vs_ref":rel(T,Tref)}
                rows.append(row)
                print(f"  r={r:.4f} rule=({nr},{nphi}) N={sc:3d} "
                      f"Eerr={row['E_rel_vs_ref']:.3e} "
                      f"Berr={row['B_rel_vs_ref']:.3e} "
                      f"Terr={row['T_rel_vs_ref']:.3e}")
            results[name][str(r)]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.15","mode":a.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: area-weighted finite-wire convergence; no GR/CTC conclusion.")

if __name__=="__main__": main()
