# AMR Brill-Lindquist Stabilitesi Raporu / AMR Brill-Lindquist Stability Report

## Türkçe Özet / Turkish Summary

### Başarı: AMR v14.1 Projeksiyonu ile Brill-Lindquist Başlangıç Verileri Stabil

**Kritik yol çözüldü.** Brill-Lindquist ikili kara deliği başlangıç verileri (zaman-simetrik, K=0) kullanılarak, AMR projeksiyon pipeline'ı (v14.1) **30+ zaman adımı için tamamen stabil** hale getirildi.

### Test Konfigürasyonu
- **Coarse Grid**: N_c=33, L_c=80.0, h_c=2.424
- **Fine Grid**: N_f=25, L_f=40.0, h_f=1.600
- **Zaman Adımı**: dt=0.03125
- **Damping**: kappa=2.0
- **Projeksiyon**: Kovaryant phi band + momentum düzeltmesi
- **Band Parametreleri**: n_layers=4, mask_threshold=0.5

### Sonuçlar

| Adımlar | Zaman | H_c_max | H_f_max | H_f_band | Durum |
|