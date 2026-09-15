
import numpy as np
C=299792458.0
G=6.67430e-11

def trace_reversed_source(T):
    eta=np.diag([1.,-1.,-1.,-1.])
    tr=np.einsum("ab,nab->n",eta,T)
    return T-0.5*eta[None,:,:]*tr[:,None,None]

def hbar_instantaneous_from_source(T,distance_m):
    Tbar=np.mean(T,axis=0)
    src=trace_reversed_source(Tbar[None,:,:])[0]
    return (4*G/C**4/max(distance_m,1e-12))*src

def characteristic_h_scale(T,size_m):
    Tbar=np.mean(np.abs(T),axis=0)
    return float((4*G/C**4)*np.sum(Tbar)*size_m)

def proper_time_fractional_shift_from_h00(h00):
    return 0.5*float(h00)
