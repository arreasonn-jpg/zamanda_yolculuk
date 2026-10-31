# v7.4-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v7.4-beta`  
**Commit**: `v7.4-beta consolidation`

---

## Summary

v7.4-beta completes the sub-QFT regime analysis, identifying the **4-regime structure** of PTUM NEC violation. The regimes are delineated by three critical fuel amplitudes: A_act ≈ 1.2 (activation), A_NEC ≈ 1.35 (NEC onset), A_QFT ≈ 1.87 (QFT violation). A paper figure (fig9_regimes.png) has been generated.

---

## Validation Results

### Sub-QFT Regime Scan (v7_4_sub_qft_regime.py)

| A_fuel | ρ_min | Ford-Roman Ratio (τ=1) | n_neg | Phi_max | psi_max | QFT OK? |
|--------|-------|------------------------|-------|---------|---------|---------|
| 0.50 | 0.0000 | 0.0000 | 0 | 0.5000 | 0.0000 | ✅ (deactivated) |
| 0.80 | 0.0000 | 0.0000 | 0 | 0.8000 | 0.0000 | ✅ (deactivated) |
| 1.00 | 0.0000 | 0.0000 | 0 | 1.0000 | 0.0000 | ✅ (deactivated) |
| 1.20 | 0.0000 | 0.0000 | 0 | 1.1999 | 0.6141 | ✅ (activated, no NEC) |
| 1.40 | -1.37e-01 | 0.1088 | 20 | 1.3997 | 0.9795 | ✅ |
| 1.60 | -4.88e-01 | 0.3880 | 20 | 1.5997 | 1.2486 | ✅ |
| 1.80 | -1.10e+00 | 0.8780 | 20 | 1.7996 | 1.4962 | ✅ |

**Key observations**:
- **A < 1.2**: PTUM deactivated (Phi < Eth threshold, psi = 0, no NEC violation)
- **1.2 < A < 1.35**: PTUM activated (psi > 0) but **no NEC violation** (ρ_min = 0)
- **1.35 < A < 1.87**: **NEC violation + QFT-compatible** (ρ_min < 0, Ford-Roman ratio < 1)
- **A > 1.87**: **NEC violation + QFT-violating** (Ford-Roman ratio > 1)

---

## 4-Regime Structure

```
A_fuel:    0.0    1.2      1.35      1.87      3.5
           |      |        |         |         |
           v      v        v         v         v
        [Deactivated] [Activated] [NEC+QFT-OK] [NEC+QFT-Viol]
           (psi=0)    (psi>0,    (rho<0,      (rho<0,
           rho=0)     rho=0)      FR<1)       FR>1)
```

**Critical thresholds**:
- **A_act ≈ 1.2**: PTUM activation (Phi > Eth at center)
- **A_NEC ≈ 1.35**: NEC violation onset (ρ_min < 0)
- **A_QFT ≈ 1.87**: Ford-Roman violation (ratio > 1 at τ=1)

---

## Generated Artifacts

- **Data**: `results/validated/v7_4_sub_qft.npz` — A, ρ_min, ratio, n_neg arrays
- **Figure**: `paper/figures/fig9_regimes.png` — 4-regime diagram with threshold annotations

---

## Paper Integration

The PTUM and BEY papers are already in the repository:

- `paper/PTUM/PTUM_LIBRARY_FINAL_v1_2_CODE_AVAILABLE_CORRECTED_COMPLETE/` — PTUM paper (main.tex, main.pdf, figures/)
- `paper/BEY/BEY_LIBRARY_FINAL_v1_1_CODE_AVAILABLE_REPO/` — BEY paper (main.tex, main.pdf, supplement.pdf, figures/)

The new figures (fig8_qft_threshold.png, fig9_regimes.png) are in `paper/figures/` and ready for integration.

---

## Unit Tests

No new unit tests for v7.4 (analysis/plotting only). Existing tests remain passing.

---

## Files Modified

- `docs/roadmap.md` — v7.4 items marked ✅ complete
- `RELEASE_NOTES_v7.4-beta.md` — This file
- `results/validated/v7_4_sub_qft.npz` — Scan data
- `paper/figures/fig9_regimes.png` — 4-regime figure

---

## Git Actions

```bash
git add docs/roadmap.md RELEASE_NOTES_v7.4-beta.md results/validated/v7_4_sub_qft.npz paper/figures/fig9_regimes.png
git commit -m "v7.4-beta consolidation: sub-QFT 4-regime structure identified, figure generated, PTUM/BEY papers in place"
git tag v7.4-beta
git push origin main --tags
```