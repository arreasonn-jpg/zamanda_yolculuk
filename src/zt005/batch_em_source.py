
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential
from .validated_stress_energy import em_stress_energy

def compute_T_batch(name, points, backpack, drive, h=0.005):
    sources=build_sources(name,backpack)
    p=np.asarray(points,float)
    n=len(p)
    omega=2*np.pi*drive.frequency_hz

    blocks=[]
    for axis in range(3):
        for delta in (2*h,h,-h,-2*h):
            q=p.copy()
            q[:,axis]+=delta
            blocks.append(q)

    Aall=vector_potential(np.vstack(blocks),sources,drive.peak_current_a)

    deriv=[]
    for axis in range(3):
        k=axis*4*n
        A2=Aall[k:k+n]
        A1=Aall[k+n:k+2*n]
        Am1=Aall[k+2*n:k+3*n]
        Am2=Aall[k+3*n:k+4*n]
        deriv.append((-A2+8*A1-8*Am1+Am2)/(12*h))

    dAx,dAy,dAz=deriv
    B=np.column_stack([
        dAy[:,2]-dAz[:,1],
        dAz[:,0]-dAx[:,2],
        dAx[:,1]-dAy[:,0],
    ])

    A0=vector_potential(p,sources,drive.peak_current_a)
    E=omega*A0
    return em_stress_energy(E,B)["T"]
