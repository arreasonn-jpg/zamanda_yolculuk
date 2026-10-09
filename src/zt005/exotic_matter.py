"""Exotic-matter stress-energy sources for CTC-relevant metrics.

Provides minimal, physically-motivated negative-energy-density models
that can be fed into the linearized GR solver alongside the EM source.

References
----------
- Morris & Thorne (1988), Am. J. Phys. 56, 395 — wormhole exotic matter
- Casimir (1948) — negative energy density between plates
- Tipler (1974), Phys. Rev. D 9, 2203 — rotating cylinder
- Gödel (1949), Rev. Mod. Phys. 21, 447 — rotating universe CTC

NOTE: this module is *not* a device model. It provides computable
stress-energy profiles whose magnitude can be compared against the EM
source to quantify the gap.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict

import numpy as np

# SI constants
C = 2.99792458e8
G = 6.67430e-11
HBAR = 1.054571817e-34


@dataclass
class ExoticSource:
    """Container for a stress-energy profile.

    Attributes
    ----------
    name : str
        Identifier used in reports.
    rho : Callable[[float], float]
        Energy density T^{00} (J/m^3). Negative for exotic matter.
    p_r : Callable[[float], float]
        Radial pressure T^{rr}.
    p_t : Callable[[float], float]
        Tangential pressure T^{theta theta}=T^{phi phi}.
    params : dict
        Input parameters (for provenance / JSON dumps).
    """

    name: str
    rho: Callable[[float], float]
    p_r: Callable[[float], float]
    p_t: Callable[[float], float]
    params: Dict[str, float]

    def sample(self, r: np.ndarray) -> Dict[str, np.ndarray]:
        r = np.asarray(r, dtype=float)
        return {
            "r": r,
            "rho": np.array([self.rho(x) for x in r]),
            "p_r": np.array([self.p_r(x) for x in r]),
            "p_t": np.array([self.p_t(x) for x in r]),
        }

    def null_energy_violation(self, r: float, k: np.ndarray) -> float:
        """T_{mu nu} k^mu k^nu for a null vector k in a local orthonormal
        frame. Returns a scalar; negative => NEC violated."""
        k = np.asarray(k, dtype=float)
        # Diagonal approximation in orthonormal frame
        rho = self.rho(r)
        p_r = self.p_r(r)
        p_t = self.p_t(r)
        return rho * k[0] ** 2 + p_r * k[1] ** 2 + p_t * (k[2] ** 2 + k[3] ** 2)


# ---------------------------------------------------------------------------
# Model 1 — Casimir vacuum between parallel plates
# ---------------------------------------------------------------------------
def casimir_parallel_plates(separation_m: float = 1e-6,
                            area_m2: float = 1.0) -> ExoticSource:
    """Negative energy density between two conducting plates.

    rho = -pi^2 hbar c / (720 d^4)  (ideal, T=0)
    Pressure is anisotropic: p_z = -3 rho (i.e. positive), p_x=p_y = -rho.
    Magnitude ~ -1.3e-3 J/m^3 at d = 1 um -> tiny but real.
    """
    d = float(separation_m)
    rho0 = -np.pi ** 2 * HBAR * C / (720.0 * d ** 4)

    def rho(r):
        return rho0

    def p_r(r):
        return -3.0 * rho0  # along plate normal

    def p_t(r):
        return -rho0

    return ExoticSource(
        name="casimir_parallel_plates",
        rho=rho,
        p_r=p_r,
        p_t=p_t,
        params={"separation_m": d, "area_m2": float(area_m2),
                "rho0_J_per_m3": rho0},
    )


# ---------------------------------------------------------------------------
# Model 2 — Morris-Thorne thin-shell exotic matter
# ---------------------------------------------------------------------------
def morris_thorne_shell(b0_m: float = 1.0,
                        throat_radius_m: float = 1.0) -> ExoticSource:
    """Simplified Morris-Thorne exotic-matter shell.

    Uses the standard shape function b(r) = b0^2 / r and the
    thin-shell approximation with rho ~ -b0^2 / (8 pi G r^4) near throat.
    """
    b0 = float(b0_m)
    r0 = float(throat_radius_m)
    pref = b0 ** 2 / (8.0 * np.pi * G)

    def rho(r):
        return -pref / (r ** 4 + 1e-30)

    def p_r(r):
        # p_r = -rho * b(r)/(r - b(r)); near throat |p_r| ~ |rho|
        b = b0 ** 2 / max(r, 1e-12)
        return -rho(r) * b / max(r - b, 1e-12)

    def p_t(r):
        return 0.5 * (rho(r) - p_r(r))

    return ExoticSource(
        name="morris_thorne_shell",
        rho=rho, p_r=p_r, p_t=p_t,
        params={"b0_m": b0, "throat_radius_m": r0},
    )


# ---------------------------------------------------------------------------
# Model 3 — Tipler rotating-cylinder exterior (near-field region)
# ---------------------------------------------------------------------------
def tipler_cylinder(radius_m: float = 1.0,
                    omega_rad_s: float = 1.0,
                    length_m: float = 10.0) -> ExoticSource:
    """Tipler cylinder energy density (weak-field exterior approximation).

    rho ~ -omega^2 r^2 / (8 pi G) * f(r/R)  inside the light-cylinder
    region. Negative only inside the CTC region.
    """
    R = float(radius_m)
    w = float(omega_rad_s)

    def rho(r):
        x = r / max(R, 1e-12)
        return -(w ** 2 * R ** 2) * (x ** 2) / (8.0 * np.pi * G)

    def p_r(r):
        return rho(r)

    def p_t(r):
        return rho(r)

    return ExoticSource(
        name="tipler_cylinder",
        rho=rho, p_r=p_r, p_t=p_t,
        params={"radius_m": R, "omega_rad_s": w, "length_m": float(length_m)},
    )


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def compare_to_em_source(em_rho_max: float, exotic: ExoticSource,
                         r_probe: float = 1.0) -> Dict[str, float]:
    """Ratio of exotic to EM energy density at r_probe.

    A ratio ~ 10^{-40} means "you are 40 orders of magnitude away".
    """
    ex_rho = exotic.rho(r_probe)
    ratio = ex_rho / max(abs(em_rho_max), 1e-300)
    return {
        "r_probe_m": r_probe,
        "em_rho_max_J_per_m3": em_rho_max,
        "exotic_rho_J_per_m3": ex_rho,
        "ratio_exotic_over_em": ratio,
        "log10_ratio": float(np.log10(abs(ratio) + 1e-300)),
    }


__all__ = [
    "ExoticSource",
    "casimir_parallel_plates",
    "morris_thorne_shell",
    "tipler_cylinder",
    "compare_to_em_source",
]
