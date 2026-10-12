import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_complex_ptum import evolve_refined_complex_ptum


def test_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params)
    assert np.allclose(out["state_c"].phi, 0.0, atol=1e-12)


def test_returns_all_keys():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=1, ptum_params=params)
    expected = ["state_c","PR_c","PI_c","PiR_c","PiI_c","Ph_c","PPh_c","ps_c","Pps_c",
                "state_f","PR_f","PI_f","PiR_f","PiI_f","Ph_f","PPh_f","ps_f","Pps_f"]
    for k in expected:
        assert k in out, f"missing {k}"


def test_small_perturbation_finite():
    N_c = 11; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    xc = np.linspace(-L_c/2, L_c/2, N_c)
    xf = np.linspace(-L_f/2, L_f/2, N_f)
    Xc,Yc,Zc = np.meshgrid(xc,xc,xc,indexing="ij")
    Xf,Yf,Zf = np.meshgrid(xf,xf,xf,indexing="ij")
    Rc2 = Xc**2+Yc**2+Zc**2; Rf2 = Xf**2+Yf**2+Zf**2
    from zt005.ptum_field import psi_equilibrium
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    PRc = 1e-5*np.exp(-Rc2/4.0); PRf = 1e-5*np.exp(-Rf2/1.0)
    Phc = 1.5*np.exp(-Rc2/8.0); Phf = 1.5*np.exp(-Rf2/2.0)
    psc = psi_equilibrium(Phc); psf = psi_equilibrium(Phf)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_complex_ptum(
        s_c, PRc, Zc, Zc, Zc, Phc, Zc, psc, Zc,
        s_f, PRf, Zf, Zf, Zf, Phf, Zf, psf, Zf,
        h_c=h_c, h_f=h_f, offset=(3,3,3),
        dt=1e-5, n_steps=2, ptum_params=params)
    for k in out:
        arr = out[k]
        if hasattr(arr, "phi"):
            assert np.all(np.isfinite(arr.phi)), f"{k} not finite"
        else:
            assert np.all(np.isfinite(arr)), f"{k} not finite"
