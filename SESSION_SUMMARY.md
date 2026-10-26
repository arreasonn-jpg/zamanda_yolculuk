# Zamanda Yolculuk - Session Summary

## Genel Durum
**Proje**: Zamanda Yolculuk - Spacetime manipulation hypothesis testing framework
**Son Güncelleme**: 26 Ekim 2026
**Test Durumu**: 779 test tamamı PASS

---

## Tamamlanan İşler (Priority Sırası)

### 1. ✅ Gauge Audit of Device Pipeline (Öncelik 1)
**Dosya**: `scripts/gauge_audit_device.py`
**Sonuç**: Tüm konfigürasyonlar (G1-G4) **FAIL**
- Scaled residual: 0.85 - 1.13 (kabul kriteri: < 0.15)
- **Kök Neden**: Magneto-quasi-static (MQS) yaklaşımı displacement current'i ihmal ediyor
- EM stress-energy tensor korunmuyor: ∂_μ T^μν ≠ 0
- Retarded linearized Einstein solver (Lorenz gauge) conserved source gerektiriyor
- Manufactured conserved source (GaussianBump) **PASS** (scaled_residual ~ 0.03)

### 2. ✅ Energy Bookkeeping Skeleton (Öncelik 2)
**Dosya**: `src/zt005/energy_budget.py`
**Sonuç**: Tüm konfigürasyonlar için hesaplandı
| Config | P_in (W) | P_loss (W) | P_stored (J) | η | h_current | Gap (orders) |
|--------|----------|------------|--------------|---|-----------|--------------|
| G1 | 84 | 84 | 7.37e-06 | 8.78e-08 | 6.09e-50 | 46.2 |
| G2 | 84 | 84 | 3.66e-06 | 4.35e-08 | 3.02e-50 | 46.5 |
| G3 | 84 | 84 | 1.85e-08 | 2.20e-10 | 1.53e-52 | 48.8 |
| G4 | 84 | 84 | 1.04e-06 | 1.24e-08 | 8.61e-51 | 47.1 |

**Çıkarım**: Hedef h=10⁻³'e ulaşmak için 46-49 orders of magnitude gap var.

### 3. ✅ BSSN AMR Projection Long Runs (Öncelik 3) - v14.1
**Dosyalar**: `src/zt005/amr/amr_projection.py`, `src/zt005/amr_long_run.py`
**Test Scriptleri**: `scripts/test_amr_*.py` (15+ test scripti)

**Bulgular**:
- **Projection çalışıyor**: Fine grid band Hamiltonian constraint violation %10-15 azalıyor
- **Coarse grid blowup**: N=17, L=24, dt=0.125 → ~17 step'te patlıyor (CFL/resolution sorunu)
- **Stabilite sınırları**:
  - Coarse grid: N_c≥17, L_c≥120, dt≤0.0625, kappa≥2.0 → 30+ step stabil
  - Fine grid: N_f≥25, L_f≥80, dt≤0.125, kappa≥2.0 → 25-30 step stabil
- **Brill-Lindquist** (time-symmetric, K=0) **Bowen-York'tan daha stabil**
- Domain size kritik: L≥120 ile N=17 stabil

### 4. ✅ BSSN ↔ Retarded_hbar Cross-Check (Öncelik 4)
**Dosya**: `src/zt005/run_zt006_3_gr_layer.py`
**Sonuç**: **16/16 check PASS**
- Metric sanity: Minkowski, Schwarzschild, Kerr, Gödel, Tipler, Morris-Thorne ✓
- Geodesic norm conservation: drift < 1e-10 ✓
- Perihelion precession: rel_err = 0.63% ✓
- CTC known-answer: Gödel, Tipler, negative controls ✓
- Conservation/Gauge: conservation_max_rel = 0.0079, lorenz_scaled_residual = 0.112 ✓

---

## Teknik Detaylar

### MQS Sınırlaması (Bilinen Sorun)
```
Magneto-quasi-static (MQS) approximation:
  B = curl A
  E = -∂A/∂t
  Ampere's law: ∇×B = μ₀J  (displacement current μ₀ε₀∂E/∂t ihmal edilmiş)
  
Sonuç: ∂_μ T^μν ≠ 0 → Lorenz gauge ∂_μ h̄^μν = 0 sağlanamıyor
```

### AMR Projection (v13.2/v13.3/v13.4 → v14.1)
```
project_phi_band_covariant(): Flat Laplacian + covariant operator
momentum_band_correct(): Atilde divergence correction (1 pass)
Parameters: n_layers=4, mask_threshold=0.5, omega=0.8, sweeps=10
```

### Stabilite Karakterizasyonu
| Config | N_c | L_c | N_f | L_f | dt | kappa | Steps | Result |
|--------|-----|-----|-----|-----|-----|-------|-------|--------|
| Bowen-York | 17 | 24 | 9 | 6 | 0.125 | 0.5 | ~17 | BLOWUP |
| Bowen-York | 17 | 48 | 17 | 16 | 0.0625 | 2.0 | 10 | Stable |
| Brill-Lindquist | 17 | 120 | 25 | 60 | 0.0625 | 2.0 | 30+ | Stable |

---

## Kalan İşler

### Hemen Yapılacaklar
1. [ ] MQS conservation limitation'ı dokümante et (known issue)
2. [ ] AMR testlerini Brill-Lindquist + larger grids ile genişlet
3. [ ] Coarse grid stability için N_c=49, L_c=80, dt=0.03125 test et
4. [ ] Fine grid outer boundary conditions araştır

### Sonraki Fazlar (Roadmap)
- **Phase 7**: Device-generated spacetime search (Phases 2/3 convergence'a bağlı)
- **ZT-007**: Full-wave Maxwell (frequency scans talep edene kadar ertelenmiş)

---

## Önemli Dosyalar

| Dosya | Açıklama |
|-------|----------|
| `scripts/gauge_audit_device.py` | Gauge audit ana script |
| `src/zt005/retarded_gr.py` | Retarded linearized Einstein solver |
| `src/zt005/gauge_checks.py` | Lorenz gauge diagnostic |
| `src/zt005/amr/amr_projection.py` | Band-local conformal projection |
| `src/zt005/amr_long_run.py` | Full AMR long run driver |
| `src/zt005/bssn_ic_binary.py` | Brill-Lindquist initial data |
| `src/zt005/energy_budget.py` | Energy bookkeeping |
| `src/zt005/run_zt006_3_gr_layer.py` | GR validation layer |
| `results/validated/gauge_audit_device.json` | Gauge audit sonuçları |
| `results/validated/energy_budget.json` | Energy budget sonuçları |
| `results/validated/zt006_3_gr_layer_results.json` | Cross-check sonuçları |

---

## Sonraki Adım Önerisi

**Öncelik**: BSSN AMR projection long runs (v14.1) ile devam et
- Brill-Lindquist initial data kullan
- Larger coarse grid (N_c=33-49, L_c=80-120)
- Smaller dt (0.03125-0.0625)
- Higher kappa (2.0-4.0)
- 30-50 step hedef

Bu configuration coarse grid blowup sorununu çözmeli ve fine grid projection'ın uzun vadeli etkisini göstermeli.