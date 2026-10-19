"""ZT-006.2.25 — retarded device-source convergence (fast regression).

Runs the smallest grid pair with one external observer and asserts the
consecutive-grid relative change stays far below the acceptance
criterion (rel < 0.05 full chains; here a loose 0.15 for a single
coarse pair). Full grids live in the stage runner, not in tests.
"""
import numpy as np

from zt005.model import ActiveCylinder, Backpack, Drive
from zt005.geometry import build_sources
from zt005.green_convergence import focused_source_nodes
from zt005.retarded_gr import retarded_hbar
from zt005.validated_stress_energy import em_stress_energy
from zt005.run_zt006_2_25 import A_WIRE_EFF, device_field_shapes


def test_retarded_device_convergence_smallest_pair():
    cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
    omega = 2.0 * np.pi * dv.frequency_hz
    obs = np.array([[0.0, 0.0, 2.0]])
    period = 1.0 / dv.frequency_hz
    times = np.array([0.0, period / 16, period / 8])

    fields = []
    for g in [(6, 12, 6), (8, 16, 8)]:  # first pair of the acceptance chain
        pts, w = focused_source_nodes(cyl.radius_m, cyl.height_m, *g)
        A1, B1 = device_field_shapes(pts, build_sources("G1", bk),
                                     dv, A_WIRE_EFF)

        def source_fn(t_ret, _p):
            tau = np.asarray(t_ret, float)[:, None]
            E = omega * A1 * np.sin(omega * tau)
            B = B1 * np.cos(omega * tau)
            return em_stress_energy(E, B)["T"]

        fields.append(retarded_hbar(source_fn, obs, times, pts, w))

    hc, hp = fields[1], fields[0]
    denom = np.linalg.norm(hc)
    if denom == 0.0:  # pragma: no cover - degenerate normalization
        raise AssertionError("retarded field norm is zero")
    rel = float(np.linalg.norm(hc - hp) / denom)
    assert np.isfinite(rel)
    assert rel < 0.15
