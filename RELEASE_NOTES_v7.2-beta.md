# v7.2-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v7.2-beta`  
**Commit**: `v7.2-beta consolidation`

---

## Summary

v7.2-beta documents the **negative result** for transient PTUM discharge: with the fiducial PTUM parameters (λ=0.05), the NEC violation **does not discharge** — Phi remains stable at ~3.0 and ρ_min stays at ~-15.8 for the full 500 steps (t=0.1). The fast-discharge variant (λ=1.0, non-fiducial) crashes due to singular matrix in the metric evolution.

---

## Validation Results

### 1. Transient Discharge (Fiducial λ=0.05, 500 steps)

| Metric | Value | Status |
|--------|-------|--------|
| **NEC violation** | ρ_min = -1.5844e+01 (persistent, slightly growing) | ❌ No discharge |
| **Phi_max** | 3.0008 → 3.0087 (stable, slightly growing) | ❌ No decay |
| **psi_max** | 2.8284 (constant) | — |
| **K_max (coarse)** | 0.2463 → 0.2738 (stable) | ✅ |
| **alpha_min** | ~0.84 (stable) | ✅ |
| **Discharge duration** | 0.0900 (entire run, ρ < -1.0) | ❌ Never recovers |

**Grid**: Coarse N=9, L=20; Fine N=9, L=10; offset=(2,2,2)  
**Time step**: dt=2e-4, 500 steps (t=0.1)  
**Runtime**: ~1180s

### 2. Ford-Roman Quantum Inequality (Full 500-step time series)

| τ | Massless Ratio | Massless Violates | Massive (μ=1) Ratio | Massive Violates |
|---|----------------|-------------------|---------------------|------------------|
| 0.01 | 0.0000 | False | 0.0000 | False |
| 0.05 | 0.0015 | False | 0.0015 | False |
| 0.10 | 0.0181 | False | 0.0184 | False |
| 0.50 | **2.9257** | **True** | **4.5714** | **True** |
| 1.00 | **23.6631** | **True** | **94.6523** | **True** |

**Interpretation**: Ford-Roman violated at τ ≥ 0.5 for both massless and massive cases. The persistent NEC violation (|ρ| ~ 15.8) is not QFT-compatible at any sampling time τ ≥ 0.5.

### 3. Fast-Discharge Variant (λ=1.0, Non-Fiducial)

**Result**: **CRASH** — Singular matrix in metric evolution at step ~200 (t=0.1)

```
RuntimeWarning: overflow encountered in exp
RuntimeWarning: invalid value encountered in multiply
numpy.linalg.LinAlgError: Singular matrix
```

**Cause**: Higher λ=1.0 causes faster Phi evolution which destabilizes the metric (gtilde becomes singular). This is expected — λ=1.0 is NOT the PTUM fiducial value (λ=0.05).

---

## Key Conclusion

**Transient discharge does not occur with fiducial PTUM parameters.** The equilibrium condition `psi = psi_equilibrium_smooth(Phi)` with one-way coupling creates a stable NEC-violating configuration that persists indefinitely. The PTUM field does not naturally decay — it requires an external trigger (not implemented).

This is a **negative result**: the "transient discharge" mechanism is not viable with the current PTUM model.

---

## Unit Tests

No new unit tests added for v7.2 (negative result documentation only). Existing v7.1 tests remain passing.

---

## Files Modified

- `docs/roadmap.md` — v7.2 items marked ✅ complete
- `RELEASE_NOTES_v7.2-beta.md` — This file
- `results/validated/v7_2_transient_discharge.npz` — Time series data (generated)

---

## Git Actions

```bash
git add docs/roadmap.md RELEASE_NOTES_v7.2-beta.md
git commit -m "v7.2-beta consolidation: transient discharge negative result documented"
git tag v7.2-beta
git push origin main --tags
```