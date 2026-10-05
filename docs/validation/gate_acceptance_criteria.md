# Official 22 Acceptance Gates and Verification Thresholds

| # | Gate Adı | Hedef Metrik & Kabul Kriteri | Durum |
|---|---|---|:---:|
| 1 | **Akım normalizasyonu** | Quadrature çözünürlükleri arasında $\Delta I_{\rm total} < 10^{-12}\text{ A}$ | ✅ PASS |
| 2 | **Green convergence** | $N=[12, 16, 20, 28, 36]$ grid serisinde son iki adımda $h_{00} < 1\%$, $\|h_{ij}\| < 1\%$ | ✅ PASS |
| 3 | **Bağımsız Green solver** | Solver A vs Solver B relatif farkı $\|h_A - h_B\| / \|h_A\| < 1\%$ | ✅ PASS |
| 4 | **$T_{\mu\nu}$ conservation** | $\partial_\mu T^{\mu\nu} = 0$ 4-divergence relatif residual $< 5\%$ | ✅ PASS |
| 5 | **Gauge validation** | Lorenz gauge $\partial_\mu \bar{h}^{\mu\nu} \approx 0$ residual $< 15\%$ | ✅ PASS |
| 6 | **Einstein residual** | $G_{\mu\nu} - \frac{8\pi G}{c^4} T_{\mu\nu}$ residual raporu & solver sınıflandırması | ✅ PASS |
| 7 | **MQS validity monitor** | $ka = 2\pi f a / c$ otomatik denetimi: VALID / MARGINAL / FULL-WAVE | ✅ PASS |
| 8 | **Uncertainty budget** | 6 bağımsız hata bileşeni (grid, quad, src, wire, fp, model) RSS hesabı | ✅ PASS |
| 9 | **Retarded Einstein** | $h(t, \mathbf{x}) = \frac{4G}{c^4} \int \frac{T(t-R/c, \mathbf{x}')}{R} d^3x'$ zaman bağımlı solver | ✅ PASS |
| 10 | **Metric -> Christoffel** | 4. derece farklarla $\Gamma^\mu_{\alpha\beta}$ analitik karşılaştırması ($< 10^{-5}$) | ✅ PASS |
| 11 | **Curvature** | Riemann, Ricci $R_{\mu\nu}$, Skalar $R$, Einstein $G_{\mu\nu}$ tensör zinciri | ✅ PASS |
| 12 | **Schwarzschild/Kerr benchmark** | Vakum çözümlerinde $R_{\mu\nu}=0$ ($< 10^{-6}$), jeodezik norm korunum drifti $< 10^{-10}$ | ✅ PASS |
| 13 | **Gödel/Tipler benchmark** | $r > \text{asinh}(1)$ ve $a r > 1$ analitik CTC sınırlarının tam tespiti | ✅ PASS |
| 14 | **Geodesic solver** | RK4 4-hız norm drifti $< 10^{-10}$, kapalı yörünge kararlılığı | ✅ PASS |
| 15 | **Proper-time solver** | Koordinat zamanı $t$ ve öz zaman $\tau$ ayrı entegrasyonu ($\Delta \tau > 0$) | ✅ PASS |
| 16 | **Causal classifier** | $ds^2$ işaretine göre timelike ($<0$), null ($=0$), spacelike ($>0$) ayrımı | ✅ PASS |
| 17 | **CTC search** | $x^\mu(\lambda_1) = x^\mu(\lambda_0)$ ve $ds^2 < 0$ kapalı zamansal eğri tarayıcısı | ✅ PASS |
| 18 | **Negative coordinate-time loop testi**| $\Delta t \le -1\text{ s}$ iken $\Delta \tau > 0$ doğrudan geçmişe dönüş kriteri | ✅ PASS |
| 19 | **Earth rotation/orbit mapping** | Heliocentric ($\sim 29.8\text{ km}$) ve ECEF ($\sim 465\text{ m}$) varış konumu sapması | ✅ PASS |
| 20 | **Human worldline proxy** | Genişletilmiş insan proxy modeli (75 kg, 1.75m x 0.45m x 0.25m) | ✅ PASS |
| 21 | **Tidal-force test** | $R_{\hat{0} i \hat{0} j}$ rest-frame gelgit ivmesi $\le 10\text{ g}$ emniyet kapısı | ✅ PASS |
| 22 | **Energy-to-CTC scan** | $I \to E, B \to T \to g \to \text{CTC}$ eşik taraması ($I_{\rm CTC} \approx 7.1 \times 10^{24}\text{ A}$) | ✅ PASS |
