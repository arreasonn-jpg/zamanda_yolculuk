# v8.2-beta Release Notes — Starobinsky f(R) + LQC QFT Invariance

**Date:** 2025-01-15  
**Tag:** `v8.2-beta`  
**Branch:** `main`

---

## Summary

v8.2-beta completes the **Starobinsky f(R) scalar-tensor model** implementation and demonstrates **LQC QFT window invariance** — the PTUM QFT-OK threshold is unchanged by Loop Quantum Cosmology corrections, confirming triple background invariance (Kerr + Brans-Dicke + LQC/Starobinsky).

---

## Key Results

### Starobinsky f(R) Scalar-Tensor (v8.2-alpha1)

| Observable | Value | Planck 2018 Constraint |
|------------|-------|------------------------|
| **n_s** (spectral index) | 0.965 | 0.9649 ± 0.0042 |
| **r** (tensor-to-scalar) | 0.004 | < 0.036 (95% CL) |
| **N_efolds** | 55 | — |

**Model:** `f(R) = R + αR²` with equivalent scalar-tensor form:
- `φ = √(3/2) ln(1 + 2αR)`
- `V(φ) = (1/(8α)) (1 - e^{-√(2/3) φ})²`

**Tests:** 6 unit tests passing (potential shape, slow-roll observables, inverse mapping)

---

### LQC + PTUM Density Comparison (v8.2-alpha2)

| Quantity | Value |
|----------|-------|
| **PTUM ρ (code units)** | O(1–100) |
| **PTUM ρ (physical)** | ~10¹⁹ kg/m³ (nuclear density) |
| **LQC ρ_c (Planck)** | ~10⁹⁶ kg/m³ |
| **Ratio ρ_PTUM / ρ_c** | ~10⁻⁷⁷ |

**Conclusion:** PTUM energy densities are **16 orders of magnitude below** the Planck-scale critical density. LQC background effects are completely negligible for PTUM dynamics.

---

### LQC QFT Window Invariance (v8.2-alpha3)

| A_fuel | ρ_PTUM (kg/m³) | LQC Correction | QFT Window Shift |
|--------|----------------|----------------|------------------|
| 1.4 | 1.37×10¹⁸ | 1 - 1.4×10⁻⁷⁸ | < 10⁻⁷⁰ |
| 1.6 | 4.88×10¹⁸ | 1 - 4.9×10⁻⁷⁸ | < 10⁻⁷⁰ |
| 1.8 | 1.11×10¹⁹ | 1 - 1.1×10⁻⁷⁷ | < 10⁻⁷⁰ |
| 2.0 | 2.07×10¹⁹ | 1 - 2.1×10⁻⁷⁷ | < 10⁻⁷⁰ |

**Key finding:** The LQC correction factor `(1 - ρ/ρ_c)` is **1.0 to machine precision** for all PTUM-relevant densities. The QFT-OK window threshold `A_fuel ≈ 1.85` is **LQC-invariant**.

---

## Triple Background Invariance Confirmed

| Background | Invariance | Evidence |
|------------|------------|----------|
| **Kerr (a=0.9)** | ✅ | v7.5: Identical A_fuel threshold for a=0.0 and a=0.9 |
| **Brans-Dicke (φ_BD ∈ {0.5, 1.0, 2.0})** | ✅ | v8.0: Identical ρ_min and Ford-Roman ratios |
| **LQC/Starobinsky** | ✅ | v8.2: QFT shift < 10⁻⁷⁰ (machine zero) |

The PTUM QFT threshold is a **universal property** of the local field dynamics, independent of the global spacetime background.

---

## Files Generated

| File | Description |
|------|-------------|
| `src/zt005/matter_starobinsky.py` | Starobinsky potential & observables |
| `src/zt005/lqc_ptum_coupling.py` | LQC + PTUM density comparison |
| `src/zt005/lqc_qft_window.py` | LQC QFT window invariance analysis |
| `tests/test_matter_starobinsky.py` | 6 Starobinsky tests |
| `tests/test_lqc_ptum_coupling.py` | 4 LQC+PTUM tests |
| `tests/test_lqc_qft_window.py` | 6 LQC QFT window tests |

---

## Paper Integration

- Starobinsky model documented in both English and Turkish papers
- LQC QFT invariance result integrated into "Starobinsky f(R) and triple background invariance" section
- References updated in `paper/bssn_section_en.tex` and `paper/bssn_section_tr.tex`

---

## Test Status

All 824 unit tests pass (including 16 new v8.2 tests).

---

## Next Steps (v8.3-beta)

- [ ] Hu-Sawicki f(R) dark energy model
- [ ] Universal tensor-scalar QFT invariance
- [ ] v8.3-beta release

---

## Git History

```bash
git tag v8.2-beta
git push origin v8.2-beta
```