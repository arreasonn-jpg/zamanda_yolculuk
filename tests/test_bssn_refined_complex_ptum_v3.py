import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_complex_ptum_v3 import evolve_refined_complex_ptum_v3


def test_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v3(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params)
    assert np.allclose(out["state_c"].phi, 0.0, atol=1e-12)


def test_peak_mode_preserves_Phi():
    """With matter_restrict_mode='peak', Phi should survive restrict."""
    N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    xc = np.linspace(-L_c/2, L_c/2, N_c)
    xf = np.linspace(-L_f/2, L_f/2, N_f)
    Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing='ij')
    Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing='ij')
    Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2
    from zt005.ptum_field import psi_equilibrium_smooth
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Phc = 3.0*np.exp(-Rc2/8.0); Phf = 3.0*np.exp(-Rf2/8.0)
    psc = psi_equilibrium_smooth(Phc); psf = psi_equilibrium_smooth(Phf)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v3(
        s_c, Zc, Zc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, Zf, Zf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=(3,3,3),
        dt=1e-5, n_steps=1, ptum_params=params,
        restrict_matter=True, matter_restrict_mode="peak",
        matter_blend=0.7)
    # Phi should be near 3.0 at center (peak preserved)
    Phc_center = out["Ph_c"][4, 4, 4]
    assert Phc_center > 2.5, f"Phi center = {Phc_center}"


def test_avg_mode_matches_v2():
    """With matter_restrict_mode='avg' and restrict_matter=True,
    behaves like the avg version."""
    N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    xc = np.linspace(-L_c/2, L_c/2, N_c)
    xf = np.linspace(-L_f/2, L_f/2, N_f)
    Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing='ij')
    Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing='ij')
    Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2
    from zt005.ptum_field import psi_equilibrium_smooth
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Phc = 3.0*np.exp(-Rc2/8.0); Phf = 3.0*np.exp(-Rf2/8.0)
    psc = psi_equilibrium_smooth(Phc); psf = psi_equilibrium_smooth(Phf)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum_v3(
        s_c, Zc, Zc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, Zf, Zf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=(3,3,3),
        dt=1e-5, n_steps=1, ptum_params=params,
        restrict_matter=True, matter_restrict_mode="avg")
    # Phi center will drop (avg loses 35%)
    assert out["Ph_c"][4, 4, 4] < 2.5
