import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_refined_sc_ptum import evolve_refined_sc_ptum
from zt005.sc_condensate import sc_condensate_potential
from zt005.kerr_sc_ptum_ic import kerr_sc_ptum_state


def test_minkowski_zero_stable():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}
    out = evolve_refined_sc_ptum(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-5, n_steps=2, V_sc=V_sc, ptum_params=params)
    assert np.allclose(out["state_c"].phi, 0.0, atol=1e-12)
    assert np.allclose(out["Phi_c"], 0.0, atol=1e-12)


def test_returns_all_keys():
    N_c = 11; N_f = 5
    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}
    out = evolve_refined_sc_ptum(
        s_c, Zc, Zc, Zc, Zc, Zc, Zc,
        s_f, Zf, Zf, Zf, Zf, Zf, Zf,
        h_c=0.5, h_f=0.25, offset=(3, 3, 3),
        dt=1e-5, n_steps=1, V_sc=V_sc, ptum_params=params)
    expected = ["state_c", "phi_sc_c", "Pi_sc_c", "Phi_c", "Pi_Phi_c",
                "psi_c", "Pi_psi_c", "state_f", "phi_sc_f", "Pi_sc_f",
                "Phi_f", "Pi_Phi_f", "psi_f", "Pi_psi_f"]
    for k in expected:
        assert k in out, f"missing {k}"


def test_small_perturbation_finite():
    N_c = 11; L_c = 20.0; h_c = L_c/(N_c-1)
    N_f = 9; L_f = 10.0; h_f = L_f/(N_f-1)
    xc = np.linspace(-L_c/2, L_c/2, N_c)
    xf = np.linspace(-L_f/2, L_f/2, N_f)
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij")
    Rc2 = Xc**2 + Yc**2 + Zc**2
    Rf2 = Xf**2 + Yf**2 + Zf**2

    s_c = minkowski_state(N_c)
    s_f = minkowski_state(N_f)
    psc_c = 1e-5 * np.exp(-Rc2/4.0)
    psc_f = 1e-5 * np.exp(-Rf2/1.0)
    Ph_c = 0.3 * np.exp(-Rc2/8.0)
    Ph_f = 0.3 * np.exp(-Rf2/2.0)
    Zc = np.zeros((N_c, N_c, N_c))
    Zf = np.zeros((N_f, N_f, N_f))
    V_sc = sc_condensate_potential()
    params = {"V0": 1.0, "lam": 0.05, "m2": 1.0, "alpha_p": 1.0,
              "Eth": 1.0, "lam_psi": 1.0}

    out = evolve_refined_sc_ptum(
        s_c, psc_c, Zc, Ph_c, Zc, Zc, Zc,
        s_f, psc_f, Zf, Ph_f, Zf, Zf, Zf,
        h_c=h_c, h_f=h_f, offset=(3, 3, 3),
        dt=1e-5, n_steps=2, V_sc=V_sc, ptum_params=params)
    for k in out:
        arr = out[k]
        if hasattr(arr, "phi"):
            assert np.all(np.isfinite(arr.phi)), f"{k} not finite"
        else:
            assert np.all(np.isfinite(arr)), f"{k} not finite"
