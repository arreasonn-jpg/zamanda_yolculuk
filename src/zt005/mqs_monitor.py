"""Magnetoquasistatic (MQS) Regime Validity Monitor.

Addresses Section 5:
Automatically calculates the dimensionless parameter:
    ka = 2 * pi * f * a / c = omega * a / c

Evaluates:
- MQS VALID          : ka < 0.1      (quasistatic assumption well satisfied)
- MQS MARGINAL       : 0.1 <= ka <= 0.5 (phase shift across apparatus non-negligible)
- FULL-WAVE REQUIRED : ka > 0.5      (retardation and radiation dominate)

Also calculates free-space wavelength lambda, device transit time vs cycle period,
and copper conductor skin depth.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0
MU0 = 4.0e-7 * np.pi
SIGMA_CU = 5.8e7  # Copper conductivity [S/m]


def evaluate_mqs_validity(frequency_hz: float, characteristic_dimension_m: float) -> dict:
    """Evaluate MQS validity and report regime diagnostics."""
    f = float(frequency_hz)
    a = float(characteristic_dimension_m)

    if f <= 0:
        raise ValueError(f"frequency must be positive, got {f}")
    if a <= 0:
        raise ValueError(f"characteristic dimension must be positive, got {a}")

    omega = 2.0 * np.pi * f
    k = omega / C
    ka = k * a
    wavelength_m = C / f
    transit_time_s = a / C
    period_s = 1.0 / f
    transit_to_period_ratio = transit_time_s / period_s

    # Skin depth in copper conductor: delta = sqrt(2 / (omega * mu0 * sigma))
    skin_depth_m = np.sqrt(2.0 / (omega * MU0 * SIGMA_CU))

    if ka < 0.1:
        status = "MQS VALID"
        explanation = "Apparatus is much smaller than EM wavelength (ka < 0.1); quasistatic approximation is strictly defensible."
    elif ka <= 0.5:
        status = "MQS MARGINAL"
        explanation = "Phase variations across apparatus are moderate (0.1 <= ka <= 0.5); full-wave treatment recommended."
    else:
        status = "FULL-WAVE REQUIRED"
        explanation = "Apparatus size approaches or exceeds EM wavelength (ka > 0.5); full-wave Maxwell solution required."

    return {
        "status": status,
        "ka": float(ka),
        "frequency_hz": f,
        "dimension_a_m": a,
        "wavelength_m": float(wavelength_m),
        "transit_time_s": float(transit_time_s),
        "period_s": float(period_s),
        "transit_to_period_ratio": float(transit_to_period_ratio),
        "skin_depth_cu_mm": float(skin_depth_m * 1e3),
        "explanation": explanation,
        "is_mqs_valid": bool(status == "MQS VALID"),
    }
