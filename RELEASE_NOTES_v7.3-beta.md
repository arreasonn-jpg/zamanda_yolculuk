# v7.3-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v7.3-beta`  
**Commit**: `v7.3-beta consolidation`

---

## Summary

v7.3-beta completes the QFT threshold scan for PTUM NEC violation. The Ford-Roman quantum inequality (massless, N=2, τ=1.0) is violated when |ρ_min| > 1.27, corresponding to A_fuel > 1.85. The scaling is linear: **ratio = 0.795 × |ρ_min|**. A paper figure (fig8_qft_threshold.png) has been generated.

---

## Validation Results

### 1. Coarse QFT Threshold Scan (v7_3_qft_threshold_scan.py)

| A_fuel | ρ_min | Ford-Roman Ratio (τ=1) | Violates |
|--------|-------|------------------------|----------|
| 1.50 | -2.84e-01 | 0.2258 | False |
| 2.00 | -2.07e+00 | 1.6442 | **True** |
| 2.50 | -6.62e+00 | 5.2639 | **True** |
| 3.00 | -1.56e+01 | 12.4308 | **True** |
| 3.50 | -3.12e+01 | 24.7885 | **True** |

**Threshold**: 1.5 < A_threshold < 2.0

### 2. Fine QFT Threshold Scan + Figure (v7_3_fine_scan_and_plot.py)

| A_fuel | ρ_min | Ford-Roman Ratio (τ=1) | Violates |
|--------|-------|------------------------|----------|
| 1.40 | -1.37e-01 | 0.1088 | False |
| 1.50 | -2.84e-01 | 0.2258 | False |
| 1.60 | -4.88e-01 | 0.3880 | False |
| 1.70 | -7.58e-01 | 0.6027 | False |
| 1.80 | -1.10e+00 | 0.8780 | False |
| **1.90** | **-1.54e+00** | **1.2222** | **True** |
| 2.00 | -2.07e+00 | 1.6442 | True |
| 2.20 | -3.47e+00 | 2.7592 | True |
| 2.50 | -6.62e+00 | 5.2639 | True |
| 3.00 | -1.56e+01 | 12.4308 | True |
| 3.50 | -3.12e+01 | 24.7885 | True |

**Precise threshold**: **A_fuel ≈ 1.85** (|ρ_min| ≈ 1.27)

### 3. Linear Scaling Law

```
Ford-Roman ratio (τ=1) = 0.7949 × |ρ_min| - 0.0002
```

R² ≈ 1.0 (perfect linear fit across all 11 data points)

**QFT-compatibility condition**: |ρ_min| < 1.27 (ratio < 1.0)

---

## Generated Artifacts

- **Data**: `results/validated/v7_3_qft_scan.npz` — A, ρ_min, ratio arrays
- **Figure**: `paper/figures/fig8_qft_threshold.png` — Two-panel figure:
  - (a) QFT violation vs NEC magnitude (|ρ_min|)
  - (b) QFT violation vs fuel amplitude (A_fuel)

---

## Key Technical Details

**Configuration**: Single-grid complex+PTUM (N=9, L=20, dt=5e-4, 100 steps)  
**PTUM params**: λ=0.05 (fiducial), m²=1.0, α_p=1.0, Eth=1.0, lam_psi=1.0  
**Ford-Roman**: Massless, N_field_components=2 (Phi + psi), τ=1.0

**Why N=2?**: PTUM has two real scalar degrees of freedom (Phi fuel field + psi proxy field). The complex scalar (Psi_R, Psi_I) is perturbative (1e-4) and doesn't contribute significantly to NEC violation.

---

## Unit Tests

No new unit tests for v7.3 (analysis/plotting only). Existing tests remain passing.

---

## Files Modified

- `docs/roadmap.md` — v7.3 items marked ✅ complete
- `RELEASE_NOTES_v7.3-beta.md` — This file
- `results/validated/v7_3_qft_scan.npz` — Scan data
- `paper/figures/fig8_qft_threshold.png` — Paper figure

---

## Git Actions

```bash
git add docs/roadmap.md RELEASE_NOTES_v7.3-beta.md results/validated/v7_3_qft_scan.npz paper/figures/fig8_qft_threshold.png
git commit -m "v7.3-beta consolidation: QFT threshold scan complete, linear scaling ratio=0.795*|rho|, figure generated"
git tag v7.3-beta
git push origin main --tags
```