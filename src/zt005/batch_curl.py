
import numpy as np
from .pointwise_curl import curl_point_4
from .maxwell import vector_potential

def curlA_batch(points, sources, current, h):
    p=np.asarray(points,float)
    displacements=[]
    for axis in range(3):
        e=np.zeros(3); e[axis]=h
        displacements.extend([p+2*e,p+e,p-e,p-2*e])
    stacked=np.concatenate(displacements,axis=0)
    A=vector_potential(stacked,sources,current)
    n=len(p)
    vals=np.array([
        -A[:n] + 8*A[n:2*n] - 8*A[2*n:3*n] + A[3*n:4*n],
        -A[4*n:5*n] + 8*A[5*n:6*n] - 8*A[6*n:7*n] + A[7*n:8*n],
        -A[8*n:9*n] + 8*A[9*n:10*n] - 8*A[10*n:11*n] + A[11*n:12*n],
    ])
    vals=vals/(12*h)
    d_dx,d_dy,d_dz=vals
    curl=np.empty_like(d_dx)
    curl[:,0]=d_dy[:,2]-d_dz[:,1]
    curl[:,1]=d_dz[:,0]-d_dx[:,2]
    curl[:,2]=d_dx[:,1]-d_dy[:,0]
    return curl

def vector_potential_batch(points,sources,current):
    return vector_potential(np.asarray(points,float),sources,current)
