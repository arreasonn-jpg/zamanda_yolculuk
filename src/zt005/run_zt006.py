
import json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy
from .linear_gr import hbar_instantaneous_from_source,characteristic_h_scale,proper_time_fractional_shift_from_h00

OUT=Path("zt006_results.json")

def sample_points(cyl,nr=5,nz=6,nphi=12,seed=0):
    rng=np.random.default_rng(seed)
    pts=[]
    for y in np.linspace(-.85*cyl.height_m/2,.85*cyl.height_m/2,nz):
        for r in np.linspace(0,.80*cyl.radius_m,nr):
            phis=[0.] if r==0 else np.linspace(0,2*np.pi,nphi,endpoint=False)
            for p in phis: pts.append([r*np.cos(p),y,r*np.sin(p)])
    extra=[]
    while len(extra)<150:
        q=rng.uniform([-.8*cyl.radius_m,-.85*cyl.height_m/2,-.8*cyl.radius_m],
                      [.8*cyl.radius_m,.85*cyl.height_m/2,.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (.8*cyl.radius_m)**2: extra.append(q)
    return np.vstack([np.asarray(pts),np.asarray(extra)])

def main():
    cyl,b,drive=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl); h=.005; omega=2*np.pi*drive.frequency_hz
    print("="*84)
    print("ZAMANDA YOLCULUK — ZT-006")
    print("LINEARIZED EINSTEIN FIELD / GRAVITATIONAL-SOURCE DIAGNOSTICS")
    print("="*84)
    print(f"Validation points: {len(pts)}")
    print("EM-only source; Earth/device background intentionally omitted.")
    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        def A_fn(p):
            return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]
        A=np.vstack([A_fn(p) for p in pts])
        B=np.vstack([curl_point_4(A_fn,p,h) for p in pts])
        E=omega*A
        se=em_stress_energy(E,B)
        T=se["T"]
        hbar=hbar_instantaneous_from_source(T,cyl.radius_m)
        h00=float(hbar[0,0]); h0i=np.asarray(hbar[0,1:],float)
        results[name]={
            "T00_mean_J_m3":float(T[:,0,0].mean()),
            "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean()),
            "Tij_mean_abs":float(np.abs(T[:,1:,1:]).mean()),
            "hbar00_local_scale":h00,
            "hbar0i_local_scale":h0i.tolist(),
            "hbar_characteristic_norm":characteristic_h_scale(T,cyl.radius_m),
            "fractional_proper_time_shift_proxy":proper_time_fractional_shift_from_h00(h00)
        }
        print(f"{name} T00={results[name]['T00_mean_J_m3']:.6e} "
              f"T0i={results[name]['T0i_mean_norm']:.6e} "
              f"h00~{h00:.6e} |h0i|~{np.linalg.norm(h0i):.6e} "
              f"dtau/t~{results[name]['fractional_proper_time_shift_proxy']:.6e}")
    OUT.write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: weak-field source diagnostics only; no CTC conclusion.")

if __name__=="__main__": main()
