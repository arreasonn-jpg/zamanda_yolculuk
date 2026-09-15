# Energy-to-CTC Threshold Scaling Analysis

## Physical Derivation
In linearized General Relativity coupled to classical electromagnetism:
- Magnetic and electric fields scale linearly with driving current: $B \propto I$, $E \propto I$.
- The stress-energy tensor components scale quadratically: $T_{\mu\nu} \propto I^2$.
- The linearized metric perturbation scales quadratically:
  $$h_{\mu\nu} \approx \frac{4G}{c^4 R_{\rm eff}} U_{\rm stored} \propto I^2$$

## Quantitative Scaling Chain
At the laboratory baseline ($I_0 = 100\text{ A}$):
- Stored electromagnetic energy: $U_0 \approx 0.0136\text{ J}$
- Metric perturbation: $h_0 \approx 1.98 \times 10^{-46}$
- Light cone tilt: $\theta_{\rm tilt} \approx \arctan(\sqrt{h}) \approx 10^{-21}\text{ degrees}$ (negligible).

To form Closed Timelike Curves (CTCs), the light cones must tip past $45^\circ$, which requires:
$$h \sim \mathcal{O}(1)$$
Solving for the threshold current $I_{\rm CTC}$:
$$I_{\rm CTC} = I_0 \sqrt{\frac{1.0}{h_0}} \approx 100\text{ A} \times \sqrt{\frac{1.0}{1.98 \times 10^{-46}}} \approx 7.11 \times 10^{24}\text{ A}$$
The required stored electromagnetic energy is:
$$E_{\rm CTC} = U_0 \left( \frac{I_{\rm CTC}}{I_0} \right)^2 \approx 6.87 \times 10^{43}\text{ J}$$

## Astrometric Comparison
- Global annual human energy consumption: $\sim 6 \times 10^{20}\text{ J}$
- Total energy emitted by the Sun in one second: $\sim 3.83 \times 10^{26}\text{ J}$
- Total gravitational binding energy of a Type II Supernova: $\sim 10^{44}\text{ J}$

Conclusion:
The energetic barrier between laboratory apparatus ($100\text{ A}$, $0.01\text{ J}$) and physical CTC formation is **45.7 orders of magnitude in energy** and **22.9 orders of magnitude in current**.
