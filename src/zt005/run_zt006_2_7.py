
import json,argparse
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy
from .stress_energy_conditioning import componentwise_errors,field_error_components

OUT=Path("zt006_2_7_results.json")

def sample_points(cyl,n=300,seed=321):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def fields_and_T(name,pts,b,d,h):
    sources=build_sources(name,b)
    omega=2*np.pi*d.frequency_hz
    def A_fn(p):
        return vector_potential(np.atleast_2d(p),sources,d.peak_current_a)[0]
    A=np.vstack([A_fn(p) for p in pts])
    B=np.vstack([curl_point_4(A_fn,p,h) for p in pts])
    E=omega*A
    T=em_stress_energy(E,B)["T"]
    return E,B,T

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,120 if args.mode=="quick" else 300)
    hs=[0.02,0.01,0.005,0.0025]
    print("="*98)
    print("ZAMANDA YOLCULUK — ZT-006.2.7")
    print("STRESS-ENERGY / FIELD COMPONENT CONDITIONING STUDY")
    print("="*98)
    print(f"Validation points: {len(pts)}")
    print("This stage separates genuine field error from full-T_mu_nu cancellation.")

    results={}
    for name in ["G1","G2","G3","G4"]:
        Eref,Bref,Tref=fields_and_T(name,pts,b,d,hs[-1])
        rows=[]
        for h in hs:
            E,B,T=fields_and_T(name,pts,b,d,h)
            rows.append({
                "h_m":h,
                "field_errors":field_error_components(E,B,Eref,Bref),
                "T_component_errors":componentwise_errors(T,Tref),
                "T00_mean":float(np.mean(T[:,0,0])),
                "T0i_mean_norm":float(np.mean(np.linalg.norm(T[:,0,1:],axis=1))),
            })
            ce=rows[-1]["T_component_errors"]
            fe=rows[-1]["field_errors"]
            print(f"{name} h={h:.4f} "
                  f"Eerr={fe['E']['relative_l2']:.3e} "
                  f"Berr={fe['B']['relative_l2']:.3e} "
                  f"T00err={ce['T00']['relative_rms']:.3e} "
                  f"T0ierr={ce['T0i']['relative_rms']:.3e} "
                  f"Tijerr={ce['Tij']['relative_rms']:.3e}")
        results[name]=rows

    OUT.write_text(json.dumps({"stage":"ZT-006.2.7","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: component-conditioning diagnostics; no GR/CTC conclusion.")

if __name__=="__main__": main()
