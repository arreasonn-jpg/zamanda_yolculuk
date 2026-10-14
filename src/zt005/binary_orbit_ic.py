"""Binary orbit initial conditions (v10.0-alpha1).

Two scalar Gaussian density peaks in circular Kepler orbit around
common center of mass. Each peak has the appropriate tangential
velocity for orbital balance.

Physics:
    v_orbit = sqrt(G_eff * M_total / (2 * a))
    omega_orb = v_orbit / (a/2)  (for equal masses)

For GW: expected f_GW = 2 * f_orbit (quadrupole, equal masses).
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11


def circular_orbit_velocity(M_total, separation):
    """Orbital velocity for two equal masses at separation a."""
    a = separation
    return np.sqrt(G * M_total / (2.0 * a))


def binary_orbit_positions(separation, omega, t=0.0):
    """Positions of the two bodies at time t."""
    r = separation / 2.0
    theta = omega * t
    return {
        "body1": (r * np.cos(theta), r * np.sin(theta), 0.0),
        "body2": (-r * np.cos(theta), -r * np.sin(theta), 0.0),
    }


def gaussian_binary_ic(N, L, A_amp, sigma, separation, omega, t=0.0,
                       vel_boost=True):
    """Two Gaussian bumps in circular orbit.

    Returns dict with 'phi_sc', 'Pi_sc', 'metadata'.
    Pi_sc is set to give each bump the correct tangential velocity.
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

    r = separation / 2.0
    th = omega * t
    c1, s1 = np.cos(th), np.sin(th)
    # Body 1 at (r*c1, r*s1, 0), Body 2 opposite
    x1, y1 = r * c1, r * s1
    x2, y2 = -r * c1, -r * s1

    R2_1 = (X - x1)**2 + (Y - y1)**2 + Z**2
    R2_2 = (X - x2)**2 + (Y - y2)**2 + Z**2

    phi_sc = A_amp * (np.exp(-R2_1/(2*sigma**2)) + np.exp(-R2_2/(2*sigma**2)))

    # Pi = dphi/dt = -v_x dphi/dx (for tangential motion)
    if vel_boost:
        # Body 1 velocity: -omega*r*sin(th), omega*r*cos(th), 0
        vx1, vy1 = -omega * r * s1, omega * r * c1
        vx2, vy2 = omega * r * s1, -omega * r * c1
        # dphi/dx = -(X-x1)/sigma^2 * bump
        dphi1_dx = -(X - x1) / sigma**2 * np.exp(-R2_1/(2*sigma**2))
        dphi1_dy = -(Y - y1) / sigma**2 * np.exp(-R2_1/(2*sigma**2))
        dphi2_dx = -(X - x2) / sigma**2 * np.exp(-R2_2/(2*sigma**2))
        dphi2_dy = -(Y - y2) / sigma**2 * np.exp(-R2_2/(2*sigma**2))
        Pi_sc = A_amp * (vx1*dphi1_dx + vy1*dphi1_dy
                          + vx2*dphi2_dx + vy2*dphi2_dy)
    else:
        Pi_sc = np.zeros((N, N, N))

    return {
        "phi_sc": phi_sc, "Pi_sc": Pi_sc,
        "body1_pos": (x1, y1), "body2_pos": (x2, y2),
        "omega": omega, "separation": separation,
    }


def expected_gw_frequency(omega_orbit):
    """For equal-mass binary in quadrupole approximation:
    f_GW = 2 * f_orbit = 2 * omega / (2 pi) = omega / pi."""
    return omega_orbit / np.pi


__all__ = [
    "circular_orbit_velocity", "binary_orbit_positions",
    "gaussian_binary_ic", "expected_gw_frequency",
]
