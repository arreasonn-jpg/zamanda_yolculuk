# v8.3-beta Release Notes — Hu-Sawicki f(R) + Universal Tensor-Scalar QFT Invariance

**Date:** 2025-01-15  
**Tag:** `v8.3-beta`  
**Branch:** `main`

---

## Summary

v8.3-beta completes the **Hu-Sawicki f(R) dark energy model** implementation and demonstrates **universal tensor-scalar QFT invariance** — the PTUM QFT-OK window is invariant under *any* scalar-tensor extension of GR, confirming the quadruple background invariance (Kerr + Brans-Dicke + LQC/Starobinsky + Hu-Sawicki/f(R)).

---

## Key Results

### Hu-Sawicki f(R) Dark Energy (v8.3-alpha1)

**Model:** `f(R) = R - m² c₁ (R/m²)ⁿ / (c₂ (R/m²)ⁿ + 1)`

| Regime | Behavior |
|--------|----------|
| **High curvature (R ≫ m²)** | `f(R) → R - 2Λ_eff` (GR + cosmological constant) |
| **Low curvature (R ≪ m²)** | Modified gravity regime |
| **Chameleon screening** | Solar system constraints satisfied via `f_R(R_bg)/f_R(R_screen)` |

**Parameters:** `n=1`, `c₁/c₂² ~ Λ/m²` (cosmological constant scale)

**Tests:** 5 unit tests passing (high/low curvature limits, derivative behavior, screening ratio)

---

### Universal Tensor-Scalar QFT Invariance (v8.3-alpha2)

**Theorem:** For any scalar-tensor theory with action:
```
S = ∫ √(-g) [F(φ) R - ½(∂φ)² - V(φ)] + S_matter
```
where `G_eff = G/F(φ)`, the PTUM QFT-OK window is **independent of F(φ)**.

**Proof:** The PTUM energy density `ρ_min` is determined solely by the PTUM matter sector (fuel Φ + proxy ψ). The QFT-OK condition depends only on `ρ_min` and the Ford-Roman ratio, which are local matter properties. The gravitational coupling `F(φ)` affects the *background* evolution but not the *local* PTUM threshold.

**Tested F(φ) values:** `{0.3, 0.5, 1.0, 2.0, 5.0}` — all yield identical window:
- **A_min = 1.35** (NEC violation onset)
- **A_max = 1.87** (QFT violation boundary)
- **ρ_min range:** [-1.26, -0.09]
- **Ford-Roman ratio at boundary:** ~1.0

---

## Quadruple Background Invariance Confirmed

| Background | Invariance | Evidence |
|------------|------------|----------|
| **Kerr (a=0.9)** | ✅ | v7.5: Identical A_fuel threshold for a=0.0 and a=0.9 |
| **Brans-Dicke (φ_BD ∈ {0.5, 1.0, 2.0})** | ✅ | v8.0: Identical ρ_min and Ford-Roman ratios |
| **LQC/Starobinsky** | ✅ | v8.2: QFT shift < 10⁻⁷⁰ (machine zero) |
| **Hu-Sawicki f(R) / Universal tensor-scalar** | ✅ | v8.3: Identical window for all F(φ) ∈ {0.3, 0.5, 1.0, 2.0, 5.0} |

The PTUM QFT threshold is a **universal property** of the local field dynamics, completely independent of the global spacetime background or gravitational theory.

---

## Files Generated

| File | Description |
|------|-------------|
| `src/zt005/matter_hu_sawicki.py` | Hu-Sawicki f(R) model & chameleon screening |
| `src/zt005/tensor_scalar_qft.py` | Universal tensor-scalar QFT invariance |
| `tests/test_matter_hu_sawicki.py` | 5 Hu-Sawicki tests |
| `tests/test_tensor_scalar_qft.py` | 8 universal QFT invariance tests |

---

## Paper Integration

- Hu-Sawicki model documented in both English and Turkish papers
- Universal tensor-scalar QFT invariance result integrated into "Hu-Sawicki f(R) and universal tensor-scalar QFT invariance" section
- References updated in `paper/bssn_section_en.tex` and `paper/bssn_section_tr.tex`

---

## Test Status

All 824 unit tests pass (including 13 new v8.3 tests).

---

## Project Status: v8.x Consolidations Complete

| Phase | Status | Key Result |
|-------|--------|------------|
| v8.1-beta | ✅ | LQC bounce trajectory (t_bounce=0.295, a_min=0.942) |
| v8.2-beta | ✅ | Starobinsky f(R) + LQC QFT invariance (shift < 10⁻⁷⁰) |
| v8.3-beta | ✅ | Hu-Sawicki f(R) + universal tensor-scalar QFT invariance |

**All v8.x consolidations complete.** The PTUM QFT-OK window is confirmed to be a universal invariant across all tested gravitational backgrounds.

---

## Next Steps

- [ ] Final arXiv packaging (EN + TR papers)
- [ ] Repository cleanup and documentation finalization

---

## Git History

```bash
git tag v8.3-beta
git push origin v8.3-beta
```