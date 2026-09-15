
import json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy, summary, angular_momentum_density

OUT=Path("zt005_4_results.json")

def sample_points(cyl,n_radial=6,n_z=7,n_phi=16,seed=0):
    rng=np.random.default_rng(seed)
    pts=[]
    for y in np.linspace(-0.9*cyl.height_m/2,0.9*cyl.height_m/2,n_z):
        for r in np.linspace(0,0.85*cyl.radius_m,n_radial):
            phis=[0.] if r==0 else np.linspace(0,2*np.pi,n_phi,endpoint=False)
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

def main():
    cyl=ActiveCylinder(); b=Backpack(); drive=Drive()
    points=sample_points(cyl)
    h=0.005
    omega=2*np.pi*drive.frequency_hz

    print("="*84)
    print("ZAMANDA YOLCULUK — ZT-005.4")
    print("VALIDATED POINTWISE E/B → STRESS-ENERGY DIAGNOSTICS")
    print("="*84)
    print(f"Validation points: {len(points)}")
    print(f"h: {h:.4f} m")
    print(f"Drive: {drive.frequency_hz:.3e} Hz, {drive.peak_current_a:.3f} A peak normalization")
    print("B is reconstructed with the validated 4th-order curl(A) stencil.")
    print("E uses the same harmonic vector-potential baseline: E_amp = omega*A0.")

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)

        def A_fn(p):
            return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]

        A=np.vstack([A_fn(p) for p in points])
        B=np.vstack([curl_point_4(A_fn,p,h) for p in points])
        # At a quarter cycle, I=I0 cos(omega t)=0 while dI/dt=-I0 omega.
        # The field phase for E and B is therefore distinct. For magnitude
        # diagnostics, evaluate E_amp = omega*A0.
        E=omega*A

        se=em_stress_energy(E,B)
        s=summary(points,E,B,se)
        results[name]=s

        j=angular_momentum_density(points,se["S"])
        J_proxy=np.linalg.norm(j,axis=1).mean()

        print(
            f"{name} Bmean={s['Bmean_T']:.6e} T "
            f"Emean={s['Emean_V_m']:.6e} V/m "
            f"u={s['u_mean_J_m3']:.6e} J/m^3 "
            f"S={s['S_mean_W_m2']:.6e} W/m^2 "
            f"Ldens_mean={J_proxy:.6e} "
            f"|KT|mean={s['Einstein_source_norm_mean_m2']:.6e} m^-2"
        )

    OUT.write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: validated pointwise EM source diagnostics; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
