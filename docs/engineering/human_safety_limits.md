# Human Extended-Body Proxy and Tidal-Force Safety Limits

## Extended-Body Human Proxy
The human body is modeled geometrically as an extended rigid body:
- Total mass: $M = 75.0\text{ kg}$
- Stature (head-to-toe extent): $H = 1.75\text{ m}$ ($\vec{\xi}_{\rm height} = (0, 0, H)$)
- Shoulder breadth: $W = 0.45\text{ m}$
- Chest depth: $D = 0.25\text{ m}$

## Geodesic Deviation and Tidal Acceleration
The relative acceleration between the traveler's head and feet is governed by the geodesic deviation equation:
$$\frac{D^2 \xi^\mu}{d\tau^2} = R^\mu{}_{\nu\alpha\beta} u^\nu u^\alpha \xi^\beta$$
Projected into the traveler's co-moving orthonormal rest frame $(e^\mu_{(0)} = u^\mu, e^\mu_{(1)}, e^\mu_{(2)}, e^\mu_{(3)})$:
$$a^i_{\rm tidal} = - c^2 R^{\hat{i}}{}_{\hat{0}\hat{j}\hat{0}} \xi^j$$

## Physiological Thresholds & Safety Gate
- **Human Physiological Acceleration Limit**:
  Maximum allowable differential acceleration across the body is $10\text{ g} \approx 98.1\text{ m/s}^2$.
- **Structural Injury / Spaghettification Limit**:
  Tidal stresses exceeding $20\text{ g}$ cause permanent skeletal and vascular rupture.
- **Maximum Tolerable Jerk**:
  $$j = \left\| \frac{d\mathbf{a}_{\rm tidal}}{d\tau} \right\| \le 50\text{ m/s}^3$$
- **Safety Gate Status**:
  * In the laboratory device field ($h \sim 10^{-46}$): tidal acceleration is $< 10^{-28}\text{ g}$ (**PASS**).
  * Near extreme gravitational curvature or black hole horizons: tidal forces trigger immediate **FAIL**.
