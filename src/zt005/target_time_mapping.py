"""Target-time kinematic mapping placeholder (PHASE 9 groundwork).

Maps a target time t_target to where the Earth (and a surface point on it)
is located.  This is currently a KINEMATIC REFERENCE-FRAME MODEL ONLY:

* circular heliocentric orbit at 1 AU (ecliptic plane, no eccentricity),
* uniform rotation with the sidereal day,
* spherical Earth, axial tilt ignored.

It does NOT include any tensor-gravity contribution of the Earth to
T_mu_nu or to the metric (see docs/assumptions.md, "Earth model").
Its purpose at this stage is purely to make the scale of the mapping
problem visible: over one day the Earth moves ~2.6e6 km along its orbit,
so any target-time specification must be tied to the Earth's position, not
to a coordinate label alone.
"""

from __future__ import annotations

import numpy as np

AU_M = 1.495978707e11
YEAR_S = 365.25 * 86400.0
SIDEREAL_DAY_S = 86164.0905
EARTH_RADIUS_M = 6.371e6
ORBIT_OMEGA = 2.0 * np.pi / YEAR_S
SPIN_OMEGA = 2.0 * np.pi / SIDEREAL_DAY_S


def earth_center_heliocentric(t_s):
    """Earth center position in heliocentric ecliptic coordinates [m]."""
    t = np.asarray(t_s, float)
    th = ORBIT_OMEGA * t
    r = np.full_like(t, AU_M)  # circular orbit: constant radius
    return np.stack([r * np.cos(th), r * np.sin(th), np.zeros_like(t)],
                    axis=-1)


def earth_rotation_angle(t_s):
    """Rotation angle of the Earth about its axis [rad]."""
    return SPIN_OMEGA * np.asarray(t_s, float)


def surface_point_earth_centered(t_s, lat_deg, lon_deg):
    """Surface point in the Earth-centered inertial frame [m].

    The frame's z-axis is the spin axis; tilt and precession are ignored
    (documented limitation).
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg) + earth_rotation_angle(t_s)
    r_cos = EARTH_RADIUS_M * np.cos(lat)
    return np.array([r_cos * np.cos(lon),
                     r_cos * np.sin(lon),
                     EARTH_RADIUS_M * np.sin(lat)])


def device_state_heliocentric(t_s, lat_deg, lon_deg):
    """Full device position: Earth center + rotating surface offset."""
    return earth_center_heliocentric(t_s) + surface_point_earth_centered(
        t_s, lat_deg, lon_deg)


def displacement_between(t1_s, t2_s, lat_deg=0.0, lon_deg=0.0):
    """How far the device's heliocentric position moves between two times."""
    d = device_state_heliocentric(t2_s, lat_deg, lon_deg) \
        - device_state_heliocentric(t1_s, lat_deg, lon_deg)
    return float(np.linalg.norm(d))
