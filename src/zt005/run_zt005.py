import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources, check_source_envelope
from .maxwell import vector_potential, biot_savart, curl_on_regular_grid, divergence_on_regular_grid
from .metrics import rel_l2, rms, max_abs

OUT=Path("results/exploratory/zt005_results")
OUT.mkdir(parents=True, exist_ok=True)

def grid(cyl,nx=25,ny=17,nz=25):
    x=np.linspace(-cyl.radius_m,cyl.radius_m,nx)
    y=np.linspace(-cyl.height_m/2,cyl.height_m/2,ny)
    z=np.linspace(-cyl.radius_m,cyl.radius_m,nz)
    X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    mask=X*X+Z*Z <= cyl.radius_m**2
    pts=np.column_stack([X.ravel(),Y.ravel(),Z.ravel()])
    return x,y,z,X,Y,Z,pts,mask

def main():
    cyl,b,drive=ActiveCylinder(),Backpack(),Drive()
    x,y,z,X,Y,Z,pts,mask=grid(cyl)
    shape=X.shape
    axes=(x,y,z)
    omega=2*np.pi*drive.frequency_hz

    print("="*72)
    print("ZAMANDA YOLCULUK — ZT-005")
    print("MAXWELL-CONSISTENT MAGNETOQUASISTATIC VALIDATION")
    print("="*72)
    print(f"Active cylinder          : D={2*cyl.radius_m:.2f} m, H={cyl.height_m:.2f} m")
    print(f"Backpack envelope        : {b.width_m:.2f} x {b.depth_m:.2f} x {b.height_m:.2f} m")
    print(f"Drive frequency          : {drive.frequency_hz:.3e} Hz")
    print(f"Peak current             : {drive.peak_current_a:.3f} A (normalization)")
    print(f"Grid                    : {shape[0]} x {shape[1]} x {shape[2]}")
    print()

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        if not check_source_envelope(sources,b):
            raise RuntimeError(f"{name} violates backpack envelope")

        A_flat=vector_potential(pts,sources,drive.peak_current_a)
        A=A_flat.reshape(*shape,3)
        B_bs=biot_savart(pts,sources,drive.peak_current_a).reshape(*shape,3)
        B_curl=curl_on_regular_grid(A,axes)

        divB=divergence_on_regular_grid(B_curl,axes)
        # Faraday test for harmonic amplitudes:
        # E(t)=omega*A0*sin(wt), B(t)=B0*cos(wt)
        # at t=pi/(4w): curl(E)=omega*curl(A0)/sqrt(2)
        # and -dB/dt=omega*B0/sqrt(2).
        E45=(omega/np.sqrt(2))*A
        curlE=curl_on_regular_grid(E45,axes)
        minus_dBdt=(omega/np.sqrt(2))*B_curl

        # Ignore one-cell boundary where finite differences have larger error.
        interior=np.zeros(shape,dtype=bool)
        interior[1:-1,1:-1,1:-1]=True
        interior &= mask

        curl_match=rel_l2(curlE,minus_dBdt,interior)
        divB_rms=rms(divB,interior)
        divB_rel=divB_rms/max(rms(B_curl,interior)/(max(abs(x[1]-x[0]),1e-30)),1e-30)
        bs_curl=rel_l2(B_curl,B_bs,interior)

        d={
            "name":name,
            "vector_potential_mean_Am":float(np.linalg.norm(A[interior],axis=-1).mean()),
            "B_biot_savart_mean_T":float(np.linalg.norm(B_bs[interior],axis=-1).mean()),
            "B_curlA_mean_T":float(np.linalg.norm(B_curl[interior],axis=-1).mean()),
            "B_curlA_vs_BiotSavart_relL2":bs_curl,
            "divB_rms_T_per_m":divB_rms,
            "divB_relative":divB_rel,
            "Faraday_curlE_vs_minus_dBdt_relL2":curl_match,
        }
        results[name]=d

        print(
            f"{name} OK "
            f"Bmean(BS)={d['B_biot_savart_mean_T']:.3e} T "
            f"Bmean(curlA)={d['B_curlA_mean_T']:.3e} T "
            f"curlA-vs-BS={bs_curl:.3e} "
            f"divB(rel)={divB_rel:.3e} "
            f"Faraday={curl_match:.3e}"
        )

        # Central y slice of |B| from curl(A)
        iy=len(y)//2
        mag=np.linalg.norm(B_curl[:,iy,:,:],axis=-1)
        plt.figure(figsize=(6,4))
        plt.imshow(mag.T,origin="lower",
                   extent=[x[0],x[-1],z[0],z[-1]],aspect="equal")
        plt.xlabel("x [m]"); plt.ylabel("z [m]")
        plt.title(f"{name} |curl A| at y=0")
        plt.colorbar(label="|B| [T]")
        plt.tight_layout()
        plt.savefig(OUT/f"{name}_B_curlA_slice.png",dpi=150)
        plt.close()

    with open("zt005_results.json","w",encoding="utf-8") as f:
        json.dump(results,f,indent=2)

    print()
    print("Results written to zt005_results.json and zt005_results/*.png")
    print("SCIENTIFIC STATUS: magnetoquasistatic Maxwell validation; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
