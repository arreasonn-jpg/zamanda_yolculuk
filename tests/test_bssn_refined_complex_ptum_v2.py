import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_complex_ptum_v2 import evolve_refined_complex_ptum_v2


def test_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v2(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params)
    assert np.allclose(out["state_c"].phi, 0.0, atol=1e-12)


def test_psi_equilibrium_preserved():
    """With equilibrium seed, psi stays near psi_equilibrium(Phi)."""
    N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    xc = np.linspace(-L_c/2, L_c/2, N_c)
    xf = np.linspace(-L_f/2, L_f/2, N_f)
    Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing="ij")
    Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing="ij")
    Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2
    from zt005.ptum_field import psi_equilibrium
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Phc = 3.0*np.exp(-Rc2/8.0); Phf = 3.0*np.exp(-Rf2/2.0)
    psc = psi_equilibrium(Phc); psf = psi_equilibrium(Phf)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v2(
        s_c, Zc, Zc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, Zf, Zf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params)
    # psi on coarse grid should NOT be near zero
    assert np.max(np.abs(out["ps_c"])) > 0.5, \
        f"ps_c max = {np.max(np.abs(out['ps_c']))}"


def test_v1_psifromeq_disabled():
    """With psi_from_equilibrium=False, behaves like v1 (psi lost)."""
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v2(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params,
        psi_from_equilibrium=False)
    assert np.allclose(out["state_c"].phi, 0.0, atol=1e-12)
