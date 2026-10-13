"""GW extraction pipeline from BSSN evolution (v9.2-alpha1).

Given BSSN state and its time derivatives, extract GW polarizations
at observer positions using the NP formalism.
"""
from __future__ import annotations
import numpy as np
from .newman_penrose import tt_polarizations, psi4_from_hddot, psi4_amplitude


def metric_perturbation_from_state(state, eta=None):
    """Extract spatial metric perturbation h_ij = gamma_ij - delta_ij."""
    if eta is None:
        eta = np.eye(3)
    gamma = np.exp(4.0 * state.phi)[..., None, None] * state.gtilde
    return gamma - eta


def extract_at_observer(state, observer_idx, k_hat, e1, e2):
    """Extract h_+, h_x at one grid point."""
    h_full = metric_perturbation_from_state(state)
    h_local = h_full[observer_idx]
    return tt_polarizations(h_local, k_hat, e1, e2)


def extract_multiple_observers(state, observers, k_hat, e1, e2):
    """Extract h_+, h_x at multiple grid points."""
    out = []
    for obs in observers:
        hp, hc = extract_at_observer(state, obs, k_hat, e1, e2)
        out.append({"obs": obs, "h_plus": hp, "h_cross": hc})
    return out


def psi4_from_history(h_plus_history, h_cross_history, dt):
    """Compute psi_4 from time series via finite differences.

    psi_4 = -1/2 (h_+ddot - i h_xddot)
    """
    hp = np.asarray(h_plus_history, float)
    hc = np.asarray(h_cross_history, float)
    if len(hp) < 3:
        return np.array([])
    hp_ddot = np.gradient(np.gradient(hp, dt), dt)
    hc_ddot = np.gradient(np.gradient(hc, dt), dt)
    return psi4_from_hddot(hp_ddot, hc_ddot)


def strain_amplitude_from_psi4_history(psi4_history, dt):
    """Estimate |h| from |psi_4| via time-domain integration."""
    psi4 = np.asarray(psi4_history, complex)
    if len(psi4) < 3:
        return np.array([])
    # Estimate omega by FFT
    freqs = np.fft.fftfreq(len(psi4), dt)
    spectrum = np.abs(np.fft.fft(psi4))
    peak_idx = np.argmax(spectrum[1:len(spectrum)//2]) + 1
    omega_peak = 2.0 * np.pi * abs(freqs[peak_idx])
    if omega_peak < 1e-30:
        return np.zeros_like(np.abs(psi4))
    return 2.0 * np.abs(psi4) / omega_peak**2


__all__ = [
    "metric_perturbation_from_state", "extract_at_observer",
    "extract_multiple_observers", "psi4_from_history",
    "strain_amplitude_from_psi4_history",
]
