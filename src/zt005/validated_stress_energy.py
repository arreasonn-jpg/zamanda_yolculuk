
import numpy as np

C = 299_792_458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7*np.pi
G = 6.67430e-11
K = 8*np.pi*G/C**4

def em_stress_energy(E,B):
    E=np.asarray(E,float); B=np.asarray(B,float)
    E2=np.sum(E*E,axis=1); B2=np.sum(B*B,axis=1)
    u=0.5*EPS0*E2 + 0.5/MU0*B2
    S=np.cross(E,B)/MU0
    sigma=(
        EPS0*(E[:,:,None]*E[:,None,:]-0.5*E2[:,None,None]*np.eye(3))
        +(1/MU0)*(B[:,:,None]*B[:,None,:]-0.5*B2[:,None,None]*np.eye(3))
    )
    T=np.zeros((len(E),4,4))
    T[:,0,0]=u
    T[:,0,1:]=S/C
    T[:,1:,0]=S/C
    T[:,1:,1:]=-sigma
    return {
        "u":u,
        "S":S,
        "sigma":sigma,
        "T":T,
        "K_T":K*T
    }

def angular_momentum_density(points,S):
    points=np.asarray(points,float)
    momentum=S/(C*C)
    return np.cross(points,momentum)

def summary(points,E,B,se):
    bmag=np.linalg.norm(B,axis=1)
    emag=np.linalg.norm(E,axis=1)
    S=np.linalg.norm(se["S"],axis=1)
    l=np.linalg.norm(angular_momentum_density(points,se["S"]),axis=1)
    kt=np.linalg.norm(se["K_T"].reshape(len(points),-1),axis=1)
    return {
        "points":len(points),
        "Bmean_T":float(bmag.mean()),
        "Bmax_T":float(bmag.max()),
        "Emean_V_m":float(emag.mean()),
        "Emax_V_m":float(emag.max()),
        "u_mean_J_m3":float(se["u"].mean()),
        "u_max_J_m3":float(se["u"].max()),
        "S_mean_W_m2":float(S.mean()),
        "S_max_W_m2":float(S.max()),
        "angular_momentum_density_mean":float(l.mean()),
        "Einstein_source_norm_mean_m2":float(kt.mean()),
        "Einstein_source_norm_max_m2":float(kt.max()),
    }
