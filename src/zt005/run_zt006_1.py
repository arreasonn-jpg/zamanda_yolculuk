
import json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .batch_curl import curlA_batch, vector_potential_batch
from .validated_stress_energy import em_stress_energy
from .spatial_gr import spatial_hbar_from_source, inverse_trace_reverse, proper_time_proxy

OUT=Path("results/exploratory/zt006_1_results.json")

def sample_uniform_cylinder(cyl,n,seed=0):
    rng=np.random.default_rng(seed)
    r=np.sqrt(rng.random(n))*cyl.radius_m
    phi=2*np.pi*rng.random(n)
    y=(rng.random(n)-0.5)*cyl.height_m
    return np.column_stack((r*np.cos(phi),y,r*np.sin(phi)))

def observation_points(cyl):
    r=cyl.radius_m
    z=.0
    return np.array([
        [0,0,0],
        [.5*r,0,0],[-.5*r,0,0],
        [0,.5*r,0],[0,-.5*r,0],
        [0,0,.5*r],[0,0,-.5*r],
        [.9*r,0,0],
    ],float)

def main():
    cyl=ActiveCylinder(); b=Backpack(); drive=Drive()
    source_n=400
    source_points=sample_uniform_cylinder(cyl,source_n,seed=42)
    obs=observation_points(cyl)
    h=.005
    omega=2*np.pi*drive.frequency_hz
    volume=np.pi*cyl.radius_m**2*cyl.height_m

    print("="*88)
    print("ZAMANDA YOLCULUK — ZT-006.1")
    print("SPATIAL LINEARIZED-GR GREEN-FUNCTION DIAGNOSTICS")
    print("="*88)
    print(f"Source samples: {source_n}")
    print(f"Observation points: {len(obs)}")
    print(f"Active volume: {volume:.6f} m^3")
    print(f"Pointwise curl h: {h:.4f} m")
    print("Source model: EM only. Earth/device background intentionally excluded.")

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        A=vector_potential_batch(source_points,sources,drive.peak_current_a)
        B=curlA_batch(source_points,sources,drive.peak_current_a,h)
        E=omega*A
        se=em_stress_energy(E,B)
        T=se["T"]

        hbar=spatial_hbar_from_source(obs,source_points,T,volume,softening_m=0.5*h)
        hphys=np.array([inverse_trace_reverse(x) for x in hbar])

        h00=hphys[:,0,0]
        h0=np.linalg.norm(hphys[:,0,1:],axis=1)
        tau=proper_time_proxy(h00)

        results[name]={
            "source_samples":source_n,
            "observation_points":obs.tolist(),
            "h00":h00.tolist(),
            "h0i_norm":h0.tolist(),
            "proper_time_fractional_proxy":tau.tolist(),
            "h00_abs_max":float(np.max(np.abs(h00))),
            "h0i_abs_max":float(np.max(h0)),
            "h00_center":float(h00[0]),
            "h0i_center":float(h0[0]),
            "tau_proxy_center":float(tau[0]),
        }

        print(f"{name} center h00={h00[0]:.6e} "
              f"|h0i|={h0[0]:.6e} "
              f"dtau/t={tau[0]:.6e} "
              f"max|h00|={np.max(np.abs(h00)):.6e} "
              f"max|h0i|={np.max(h0):.6e}")

    OUT.write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: quasi-static linearized-GR Green-function diagnostic; no full GR/CTC conclusion.")

if __name__=="__main__": main()
