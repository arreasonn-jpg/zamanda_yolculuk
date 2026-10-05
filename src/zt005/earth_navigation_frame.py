"""Earth Navigation Reference Frame & Target-Time Coordinate Mapping.

Addresses Section 14 and 15:
- Replaces static "device motion = false" assumption with physical reference frames:
    1. Heliocentric Ecliptic Frame
    2. Earth-Centered Inertial (ECI) Frame
    3. Earth-Centered Earth-Fixed (ECEF, WGS-84) Frame
    4. Local Topocentric Frame (East, North, Up)
- Models:
    * Earth orbital motion: v_orb ~ 29.784 km/s around Sun
    * Earth sidereal rotation: omega_earth = 7.2921159e-5 rad/s
    * Active region orientation tracking
- Solves arrival location Delta x_arrival for target time displacement Delta t:
    Shows that for Delta t = -1.0 s, the Earth orbital displacement is ~29.8 km
    and Earth rotational surface displacement is ~465*cos(lat) meters.
"""

from __future__ import annotations

import numpy as np

# Physical and Geodetic Constants (WGS-84 / Astronomical)
AU_M = 1.495978707e11
ORBIT_OMEGA = 2.0 * np.pi / (365.256363 * 86400.0)  # Sidereal year [rad/s]
ORBIT_SPEED_M_S = ORBIT_OMEGA * AU_M                 # ~ 29,784 m/s

EARTH_A_M = 6378137.0                               # WGS-84 semi-major axis [m]
EARTH_F = 1.0 / 298.257223563                       # Flattening
EARTH_B_M = EARTH_A_M * (1.0 - EARTH_F)             # Semi-minor axis [m]
EARTH_E2 = 1.0 - (EARTH_B_M / EARTH_A_M)**2         # First eccentricity squared
SPIN_OMEGA = 7.2921159e-5                           # Earth sidereal rotation [rad/s]


def geodetic_to_ecef(lat_deg: float, lon_deg: float, alt_m: float) -> np.ndarray:
    """Convert WGS-84 geodetic coordinates (lat, lon, alt) to ECEF Cartesian vector [m]."""
    phi = np.deg2rad(lat_deg)
    lam = np.deg2rad(lon_deg)
    h = float(alt_m)

    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    N = EARTH_A_M / np.sqrt(1.0 - EARTH_E2 * sin_phi**2)

    x = (N + h) * cos_phi * np.cos(lam)
    y = (N + h) * cos_phi * np.sin(lam)
    z = (N * (1.0 - EARTH_E2) + h) * sin_phi
    return np.array([x, y, z], dtype=float)


def ecef_to_geodetic(xyz: np.ndarray) -> tuple[float, float, float]:
    """Convert ECEF Cartesian coordinates [m] to WGS-84 geodetic (lat_deg, lon_deg, alt_m)."""
    x, y, z = xyz
    p = np.sqrt(x*x + y*y)
    lam = np.arctan2(y, x)

    # Bowring's closed-form iterative algorithm
    theta = np.arctan2(z * EARTH_A_M, p * EARTH_B_M)
    phi = np.arctan2(
        z + EARTH_E2 * (1.0 - EARTH_E2)**(-0.5) * EARTH_A_M * np.sin(theta)**3,
        p - EARTH_E2 * EARTH_A_M * np.cos(theta)**3
    )
    sin_phi = np.sin(phi)
    N = EARTH_A_M / np.sqrt(1.0 - EARTH_E2 * sin_phi**2)
    h = p / np.cos(phi) - N

    return float(np.rad2deg(phi)), float(np.rad2deg(lam)), float(h)


def heliocentric_earth_position(t_seconds: float) -> np.ndarray:
    """Earth center in heliocentric frame at coordinate time t [m]."""
    th = ORBIT_OMEGA * t_seconds
    return np.array([AU_M * np.cos(th), AU_M * np.sin(th), 0.0], dtype=float)


def heliocentric_earth_velocity(t_seconds: float) -> np.ndarray:
    """Earth orbital velocity vector at coordinate time t [m/s]."""
    th = ORBIT_OMEGA * t_seconds
    return np.array([-ORBIT_SPEED_M_S * np.sin(th), ORBIT_SPEED_M_S * np.cos(th), 0.0], dtype=float)


def compute_arrival_mapping(
    initial_lat_deg: float = 41.0082,   # Default: Istanbul / test laboratory
    initial_lon_deg: float = 28.9784,
    initial_alt_m: float = 50.0,
    target_dt_seconds: float = -1.0,    # Target jump into past (-1.0 s)
    t0_seconds: float = 0.0,
    co_moving_lock: bool = False,       # Whether apparatus has an active Earth-frame lock
) -> dict:
    """Compute physical displacement Delta x_arrival across celestial reference frames.

    Addresses Section 14 and 15:
    For Delta t = -1.0 s:
    - In heliocentric frame: Earth was 29.78 km behind its current orbital position.
    - In Earth-centered inertial frame: Earth surface has rotated backwards by omega * Delta t.
    """
    dt = float(target_dt_seconds)
    t1 = t0_seconds + dt

    # 1. Initial positions
    p_ecef_0 = geodetic_to_ecef(initial_lat_deg, initial_lon_deg, initial_alt_m)
    p_earth_helio_0 = heliocentric_earth_position(t0_seconds)
    p_device_helio_0 = p_earth_helio_0 + p_ecef_0

    # 2. Orbital shift of Earth center
    p_earth_helio_1 = heliocentric_earth_position(t1)
    v_orbit_mean = heliocentric_earth_velocity(t0_seconds)
    delta_x_orbit = p_earth_helio_1 - p_earth_helio_0
    orbit_displacement_norm = float(np.linalg.norm(delta_x_orbit))

    # 3. Earth rotation shift: angle delta_theta = omega * dt
    delta_theta = SPIN_OMEGA * dt
    rot_matrix = np.array([
        [np.cos(delta_theta), -np.sin(delta_theta), 0.0],
        [np.sin(delta_theta),  np.cos(delta_theta), 0.0],
        [0.0,                  0.0,                 1.0],
    ])
    p_eci_rotated = rot_matrix @ p_ecef_0
    delta_x_spin = p_eci_rotated - p_ecef_0
    spin_displacement_norm = float(np.linalg.norm(delta_x_spin))

    # 4. Heliocentric displacement without co-moving lock:
    # Traveler emerges at original heliocentric position while Earth is at p_earth_helio_1
    p_device_helio_arrival = p_device_helio_0 if not co_moving_lock else (p_earth_helio_1 + p_ecef_0)
    delta_x_helio_total = p_device_helio_arrival - (p_earth_helio_1 + p_ecef_0)
    total_helio_miss_m = float(np.linalg.norm(delta_x_helio_total))

    # 5. Surface arrival location if arriving in Earth coordinates:
    if co_moving_lock:
        arr_lat, arr_lon, arr_alt = initial_lat_deg, initial_lon_deg, initial_alt_m
        surface_miss_m = 0.0
    else:
        # If emerging in ECI frame at uncompensated angle:
        arr_lat, arr_lon, arr_alt = ecef_to_geodetic(p_eci_rotated)
        surface_miss_m = spin_displacement_norm

    return {
        "initial_geodetic": {"lat_deg": initial_lat_deg, "lon_deg": initial_lon_deg, "alt_m": initial_alt_m},
        "target_dt_s": dt,
        "co_moving_lock": co_moving_lock,
        "earth_orbital_speed_m_s": float(ORBIT_SPEED_M_S),
        "earth_orbital_displacement_m": orbit_displacement_norm,
        "earth_spin_displacement_m": spin_displacement_norm,
        "total_heliocentric_displacement_m": total_helio_miss_m,
        "arrival_geodetic": {"lat_deg": arr_lat, "lon_deg": arr_lon, "alt_m": arr_alt},
        "surface_displacement_m": surface_miss_m,
        "physical_interpretation": (
            f"Over {dt:.2f} s, Earth moves {orbit_displacement_norm/1e3:.2f} km along its orbit "
            f"and the surface rotates {spin_displacement_norm:.1f} m. "
            f"{'Frame lock active: traveler remains on apparatus.' if co_moving_lock else 'No frame lock: traveler arrives 29.8 km off in space!'}"
        ),
    }
