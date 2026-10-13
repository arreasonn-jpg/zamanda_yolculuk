"""Brans-Dicke observational predictions (v8.0-alpha3).

Mercury perihelion precession:
  GR: Delta_omega_GR = 6 pi G M / (a (1 - e^2) c^2)
  BD: Delta_omega_BD = (4 + 3 omega) / (6 + 4 omega) * Delta_omega_GR

For omega = 40000 (Cassini bound), BD/GR = 1 + epsilon with
epsilon ~ 3/(4 omega) ~ 1.9e-5.
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
M_SUN = 1.989e30
AU = 1.496e11


def gr_perihelion_precession(M, a, e):
    """GR perihelion precession per orbit (radians)."""
    return 6.0 * np.pi * G * M / (a * (1.0 - e**2) * C**2)


def bd_perihelion_precession(M, a, e, omega):
    """BD-corrected precession per orbit."""
    factor = (4.0 + 3.0 * omega) / (6.0 + 3.0 * omega)
    return factor * gr_perihelion_precession(M, a, e)


def bd_deviation(omega):
    """Fractional deviation from GR."""
    factor = (4.0 + 3.0 * omega) / (6.0 + 3.0 * omega)
    return factor - 1.0


def mercury_test(omega=40000.0):
    """Mercury perihelion: GR vs BD, arcsec per century."""
    a_Mercury = 0.387 * AU
    e_Mercury = 0.2056
    P_years = 0.2408  # orbital period

    gr_rad = gr_perihelion_precession(M_SUN, a_Mercury, e_Mercury)
    bd_rad = bd_perihelion_precession(M_SUN, a_Mercury, e_Mercury, omega)

    # Convert to arcsec/century
    orbit_to_arcsec = 180.0 / np.pi * 3600.0
    orbits_per_century = 100.0 / P_years

    gr_arcsec = gr_rad * orbit_to_arcsec * orbits_per_century
    bd_arcsec = bd_rad * orbit_to_arcsec * orbits_per_century

    return {
        "gr_arcsec_per_century": gr_arcsec,
        "bd_arcsec_per_century": bd_arcsec,
        "deviation_arcsec": bd_arcsec - gr_arcsec,
        "deviation_fraction": bd_rad / gr_rad - 1.0,
    }


def gw_speed_ratio(omega):
    """GW speed in BD vs GR.
    c_GW/c ~ 1 - 1/(2 + omega) approximately for scalar-tensor.
    """
    return 1.0 - 1.0 / (2.0 + omega)


__all__ = [
    "gr_perihelion_precession", "bd_perihelion_precession",
    "bd_deviation", "mercury_test", "gw_speed_ratio",
]
