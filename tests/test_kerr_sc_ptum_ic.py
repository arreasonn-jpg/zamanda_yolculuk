import numpy as np
from zt005.kerr_sc_ptum_ic import kerr_sc_ptum_state, combined_diagnostics


def test_ic_shapes():
    pack = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0)
    assert pack["state"].phi.shape == (9,9,9)
    assert pack["phi_sc"].shape == (9,9,9)
    assert pack["Phi"].shape == (9,9,9)
    assert pack["psi"].shape == (9,9,9)


def test_ic_phi_sc_small():
    pack = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0,
                               phi_amp_target=1e-4)
    assert np.max(np.abs(pack["phi_sc"])) < 1e-3


def test_ic_ptum_activated():
    """A_fuel=3.0 should activate some cells."""
    pack = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0,
                               A_fuel=3.0)
    assert np.sum(pack["psi"] > 0) > 0


def test_diagnostics_structure():
    pack = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0)
    d = combined_diagnostics(pack, h=20.0/8.0)
    for k in ("rho_total_min", "rho_sc_min", "rho_ptum_min",
              "n_neg_total", "Kmax", "alpha_min"):
        assert k in d
    assert np.isfinite(d["rho_total_min"])


def test_total_is_sum_of_parts():
    pack = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.5, rc=2.0)
    d = combined_diagnostics(pack, h=20.0/8.0)
    # PTUM dominates (A_fuel=3.0 -> rho ~ -0.2)
    assert d["rho_ptum_min"] < 0
    # Total min should be <= SC min and <= PTUM min (superposition)
    assert d["rho_total_min"] <= d["rho_sc_min"] + 1e-6
    assert d["rho_total_min"] <= d["rho_ptum_min"] + 1e-6


def test_a_dependence():
    """Kerr spin changes beta pattern (KS has radial shift even at a=0)."""
    p0 = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.0, rc=2.0)
    p9 = kerr_sc_ptum_state(N=9, L=20.0, M=1.0, a=0.9, rc=2.0)
    assert not np.allclose(p0["state"].beta, p9["state"].beta)
