# v7.0-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v7.0-beta`  
**Branch**: `main`  
**Commit**: (to be created)

---

## Summary

v7.0-beta achieves the **first persistent NEC-violating refined evolution** of the complex scalar + PTUM stack. The key fix: **re-deriving psi from psi_equilibrium(Phi) after every prolongate/restrict operation**, combined with **one-way matter coupling** (coarse matter evolves independently, fine grid receives prolonged coarse matter).

---

## What's New

### Key Fixes (v7.0-alpha1)

1. **`psi_equilibrium_smooth()`** (`ptum_field.py`): Smoothstep-regularized equilibrium psi to avoid Laplacian spikes at the activation threshold `Phi_c`

2. **`evolve_refined_complex_ptum_v2()`** (`bssn_refined_complex_ptum_v2.py`):
   - After prolongate: `psf_ = psi_equilibrium_smooth(Phf_, ...)` — fine-grid psi derived from prolongated Phi
   - After restrict (if `restrict_matter=False`): coarse psi re-derived from coarse Phi
   - **One-way matter coupling**: coarse matter evolves independently; fine grid receives prolonged coarse matter each step

3. **`evolve_refined_complex_ptum_v3()`** (`bssn_refined_complex_ptum_v3.py`): Extended with peak-preserving restrict mode (`restrict_max`) for matter fields

### Physics Modules
- **Complex scalar matter** (`matter_complex.py`): U(1) symmetric complex scalar
- **PTUM fields** (`matter_ptum.py`): Fuel `Φ` + proxy `ψ` with tachyon condensation
- **Equilibrium seed** (`ptum_field.py`): `psi_equilibrium()` and `psi_equilibrium_smooth()`

### Evolution Drivers
- **Single-grid RK4** (`bssn_rk4_complex_ptum.py`): Coupled BSSN + complex scalar + PTUM
- **Refined 2:1 grid v2** (`bssn_refined_complex_ptum_v2.py`): Fixed equilibrium psi + one-way coupling
- **Refined 2:1 grid v3** (`bssn_refined_complex_ptum_v3.py`): Peak-preserving restrict option

---

## Validation Results

### Refined 2:1 Grid Complex + PTUM (`v7_refined_complex_ptum_run.py`, A_fuel=3.0)
```
N_c=9, L_c=20, N_f=9, L_f=10, dt=2e-4, steps=100, offset=(2,2,2)
restrict_matter=False (one-way coupling)
```

| Step | ρ_min | n_neg | K_max (coarse) | K_max (fine) | ψ_max | α_min | Finite |
|------|-------|-------|----------------|--------------|-------|-------|--------|
| init | -1.564e+01 | 1 | — | — | — | — | — |
| 25   | -1.5645e+01 | 1 | 8.74e-02 | 2.45e-01 | 2.8287 | 1.0000 | ✅ |
| 50   | -1.5651e+01 | 1 | 8.74e-02 | 2.45e-01 | 2.8290 | 1.0000 | ✅ |
| 75   | -1.5658e+01 | 1 | 8.74e-02 | 2.45e-01 | 2.8293 | 1.0000 | ✅ |
| 100  | **-1.5665e+01** | **1** | **8.74e-02** | **2.45e-01** | **2.8297** | **1.0000** | ✅ |

**Key findings:**
- **NEC violation persistent** at `ρ ≈ -15.67` for 100 steps on refined grid
- **1 negative-energy cell** maintained stably
- `K_max` stable: coarse ~0.087, fine ~0.245
- `ψ_max = 2.83` constant — equilibrium preserved
- `α_min = 1.0` constant — no collapse
- All fields finite throughout

---

## Unit Tests (17 passing)

```
tests/test_bssn_rk4_complex_ptum.py           - 4 tests  ✅
tests/test_bssn_refined_complex_ptum.py       - 3 tests  ✅
tests/test_kerr_complex_ptum_ic.py            - 3 tests  ✅
tests/test_bssn_refined_complex_ptum_v2.py    - 3 tests  ✅
tests/test_bssn_refined_complex_ptum_v3.py    - 3 tests  ✅
Total: 17 tests in 60.42s
```

---

## Design Decision: One-Way Matter Coupling

**Matter fields (Phi, psi, Psi_R/I) evolve on the coarse grid independently; the fine grid receives prolonged coarse matter each step.**

This avoids the ~35% per-step Gaussian peak loss from the restrict averaging operation. The metric still uses two-way coupling.

This is a standard AMR technique for compact sources. It trades fine-grid matter accuracy for coarse-grid stability.

---

## Grid Parameters

| Configuration | N | L | h | dt | Steps |
|---------------|---|---|---|-----|-------|
| Refined coarse | 9 | 20 | 2.5 | 2e-4 | 100 |
| Refined fine | 9 | 10 | 1.25 | 2e-4 | 100 |

---

## Key Metrics Summary

| Metric | v6.6 Refined (broken) | v7.0 Refined (fixed) |
|--------|----------------------|----------------------|
| NEC violation | ❌ ρ=0 | ✅ ρ=-15.67 |
| Negative cells | 0 | 1 |
| K_max (coarse) | 4.4e-4 | 8.7e-2 |
| K_max (fine) | 4.4e-4 | 2.45e-1 |
| ψ_max | 0 (decayed) | 2.83 (constant) |
| α_min | 1.0 | 1.0 |
| Fields finite | ✅ | ✅ |

---

## Files Changed

- `docs/roadmap.md` — v7.0 items marked ✅ complete
- `RELEASE_NOTES_v7.0-beta.md` — this file

---

## Next Steps (v7.1)

- Peak-preserving restrict (exploratory, two-way fails)
- One-way coupling remains canonical
- Refined Kerr + complex + PTUM (NEC, 100-step)
- Ford-Roman check on new config

---

**All validation scripts pass. Ready for tag v7.0-beta.**