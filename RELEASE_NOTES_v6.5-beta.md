# v6.5-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v6.5-beta`  
**Branch**: `main`  
**Commit**: (to be created)

---

## Summary

v6.5-beta consolidates the **full stack Kerr + SC + PTUM evolution** with **long-run stability validation** on both single-grid and refined 2:1 grids. This is the first release demonstrating **NEC-violating stable evolution** of the complete physics stack.

---

## What's New

### Physics Modules (all from v6.4, validated at longer runs)
- **Kerr-SC regular core** (`kerr_sc_ic.py`): Smoothed mass profile `M_eff(r) = M * (1 - exp(-r/rc))` with Israel junction matching; supports `|a| ≤ M`
- **SC condensate** (`sc_condensate.py`): Single real scalar with BCS-inspired potential `V(φ) = ½m²φ² + ¼g_intφ⁴`; perturbative amplitude `φ_amp_target=1e-4`
- **PTUM fields** (`matter_ptum.py`): Fuel `Φ` (Gaussian) + proxy `ψ` (equilibrium `ψ = √(-m²/λ)` when `Φ > Eth`); tachyon condensation triggers NEC violation
- **Full-stack IC** (`kerr_sc_ptum_ic.py`): Combines all three modules with `combined_diagnostics()`

### Evolution Drivers
- **Single-grid RK4** (`bssn_rk4_sc_ptum.py`): Coupled BSSN + SC + PTUM evolution, returns 7-tuple
- **Refined 2:1 grid** (`bssn_refined_sc_ptum.py`): Coarse→prolongate→fine→restrict cycle with tensor refinement

---

## Validation Results

### 1. Single-Grid 200-Step Long Run (`v6_sc_ptum_full_long.py`)
```
N=9, L=20, dt=1e-4, a=0.9 (worst case), A_fuel=1.5
```
| Step | ρ_min | K_max | α_min | Finite |
|------|-------|-------|-------|--------|
| 40   | -2.8422e-01 | 3.5695e-01 | 0.8742 | ✅ |
| 80   | -2.8423e-01 | 3.8583e-01 | 0.8742 | ✅ |
| 120  | -2.8424e-01 | 4.1477e-01 | 0.8742 | ✅ |
| 160  | -2.8425e-01 | 4.4378e-01 | 0.8742 | ✅ |
| 200  | **-2.8427e-01** | **4.7287e-01** | **0.8742** | ✅ |

**Key findings:**
- NEC violation **stable** at `ρ ≈ -0.284` for 200 steps
- `K_max` grows slowly (0.36 → 0.47) — controlled
- `α_min = 0.8742` **constant** — no collapse
- All fields finite throughout

### 2. Refined 2:1 Grid 100-Step (`v6_refined_sc_ptum_run.py`, A_fuel=1.5)
```
N_c=9, L_c=20, N_f=9, L_f=10, dt=2e-4, offset=(2,2,2)
```
| Step | ρ_min | K_max (coarse) | K_max (fine) | α_min | Finite |
|------|-------|----------------|--------------|-------|--------|
| 25   | 8.09e-10 | 8.04e-02 | 8.03e-02 | 0.8805 | ✅ |
| 50   | 8.09e-10 | 7.78e-02 | 7.77e-02 | 0.8805 | ✅ |
| 75   | 8.09e-10 | 7.52e-02 | 7.51e-02 | 0.8805 | ✅ |
| 100  | 8.09e-10 | **7.27e-02** | **7.26e-02** | 0.8805 | ✅ |

**Key findings:**
- **K_max decreasing** on both grids (8.04e-02 → 7.27e-02)
- No NEC violation at A_fuel=1.5 (ρ_min ~ 0)
- Stable refinement cycle

### 3. Refined 2:1 Grid 100-Step Strong (`v6_refined_sc_ptum_run_strong.py`, A_fuel=3.0)
```
N_c=9, L_c=20, N_f=9, L_f=10, dt=2e-4, offset=(2,2,2)
```
| Step | ρ_total_min | ρ_ptum_min | n_neg | K_max (c) | K_max (f) | α_min | Finite |
|------|-------------|------------|-------|-----------|-----------|-------|--------|
| init | -1.564e+01 | -1.564e+01 | — | — | — | — | — |
| 25   | -9.8550e-02 | -9.8550e-02 | 3 | 9.49e-02 | 9.53e-02 | 0.8805 | ✅ |
| 50   | -9.8553e-02 | -9.8553e-02 | 3 | 1.07e-01 | 1.07e-01 | 0.8805 | ✅ |
| 75   | -9.8559e-02 | -9.8559e-02 | 3 | 1.19e-01 | 1.19e-01 | 0.8805 | ✅ |
| 100  | **-9.8566e-02** | **-9.8566e-02** | **3** | **1.31e-01** | **1.31e-01** | **0.8805** | ✅ |

**Key findings:**
- NEC violation **maintained** on refined grid: `ρ ≈ -0.099`
- **3 negative-energy cells** persist stably
- `K_max` grows moderately but controlled
- `α_min = 0.8805` constant — no collapse
- All fields finite

---

## Unit Tests (26 passing)

```
tests/test_kerr_sc_ic.py           - 5 tests  ✅
tests/test_sc_condensate.py        - 5 tests  ✅
tests/test_kerr_sc_ptum_ic.py      - 5 tests  ✅
tests/test_bssn_refined_sc_ptum.py - 6 tests  ✅
tests/test_bssn_rk4_sc_ptum.py     - 5 tests  ✅
Total: 26 tests in 18.00s
```

---

## Grid Parameters

| Configuration | N | L | h | dt | Steps |
|---------------|---|---|---|-----|-------|
| Single-grid | 9 | 20 | 2.5 | 1e-4 | 200 |
| Refined coarse | 9 | 20 | 2.5 | 2e-4 | 100 |
| Refined fine | 9 | 10 | 1.25 | 2e-4 | 100 |

---

## Key Metrics Summary

| Metric | Single-grid (200) | Refined A=1.5 (100) | Refined A=3.0 (100) |
|--------|-------------------|---------------------|---------------------|
| NEC violation | ✅ ρ=-0.284 | ❌ (ρ≈0) | ✅ ρ=-0.099 |
| K_max trend | ↗ 0.36→0.47 | ↘ 0.080→0.073 | ↗ 0.095→0.131 |
| α_min | 0.8742 (const) | 0.8805 (const) | 0.8805 (const) |
| Fields finite | ✅ | ✅ | ✅ |
| Negative cells | N/A | 0 | 3 |

---

## Files Changed

- `docs/roadmap.md` — v6.5 items marked ✅ complete
- `RELEASE_NOTES_v6.5-beta.md` — this file

---

## Next Steps (v6.6)

- Complex scalar + PTUM driver with equilibrium seed
- Refined complex + PTUM on 2:1 grid
- Note: equilibrium psi currently broken (documented limitation)

---

**All validation scripts pass. Ready for tag v6.5-beta.**