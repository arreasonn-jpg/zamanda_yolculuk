import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_kerr_complex_ptum import evolve_refined_kerr_complex_ptum


def test_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c); s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c,N_c,N_c)); Zf = np.zeros((N_f,N_f,N_f))
    params = {"V0":1.0,"lam":0.05,"m2":1.0,"alpha_p":1.0,"Eth":1.0,"lam_psi":1.0}
    out = evolve_refined_kerr_complex_ptum(
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
    out = evolve_refined_kerr_complex_ptum(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3,3,3),
        dt=1e-5, n_steps=1, ptum_params=params)
    for k in ("state_c","Ph_c","ps_c","state_f","Ph_f","ps_f"):
        assert k in out, f"missing {k}"


def test_kerr_perturbation_finite():
    from zt005.kerr_complex_ptum_ic import kerr_complex_ptum_state
    N_c = 9; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    p_c = kerr_complex_ptum_state(N=N_c, L=L_c, M=1.0, a=0.5, A_fuel=3.0,
                                    sigma_fuel=4.0)
    p_f = kerr_complex_ptum_state(N=N_f, L=L_f, M=1.0, a=0.5, A_fuel=3.0,
                                    sigma_fuel=4.0)
    params = p_c["ptum_params"]
    out = evolve_refined_kerr_complex_ptum(
        p_c["state"], p_c["PR"], p_c["PI"], p_c["PiR"], p_c["PiI"],
        p_c["Phi"], p_c["Pi_Phi"], p_c["psi"], p_c["Pi_psi"],
        p_f["state"], p_f["PR"], p_f["PI"], p_f["PiR"], p_f["PiI"],
        p_f["Phi"], p_f["Pi_Phi"], p_f["psi"], p_f["Pi_psi"],
        h_c=h_c, h_f=h_f, offset=(3,3,3),
        dt=1e-5, n_steps=1, ptum_params=params)
    for k in out:
        arr = out[k]
        if hasattr(arr, "phi"):
            assert np.all(np.isfinite(arr.phi)), f"{k} not finite"
        else:
            assert np.all(np.isfinite(arr)), f"{k} not finite"
