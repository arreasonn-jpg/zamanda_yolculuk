# Stress-Energy Tensor Physics and Energy Conditions

## 1. Electromagnetic Stress-Energy Tensor
In SI units, the symmetric stress-energy tensor of the classical electromagnetic field is:
$$T^{00} = u = \frac{1}{2} \left( \epsilon_0 E^2 + \frac{B^2}{\mu_0} \right)$$
$$T^{0i} = \frac{S^i}{c} = \frac{1}{c\mu_0} (\mathbf{E} \times \mathbf{B})^i$$
$$T^{ij} = -\sigma^{ij} = -\left[ \epsilon_0 \left( E^i E^j - \frac{1}{2}\delta^{ij} E^2 \right) + \frac{1}{\mu_0} \left( B^i B^j - \frac{1}{2}\delta^{ij} B^2 \right) \right]$$

## 2. Four-Divergence Conservation
Local conservation of four-momentum requires:
$$\partial_\mu T^{\mu\nu} = 0 \quad (\nu = 0, 1, 2, 3)$$
- $\nu = 0$: Local energy conservation (Poynting's theorem):
  $$\frac{\partial u}{\partial t} + \nabla \cdot \mathbf{S} = -\mathbf{J} \cdot \mathbf{E}$$
- $\nu = i$: Local momentum conservation (Maxwell stress):
  $$\frac{1}{c^2} \frac{\partial \mathbf{S}}{\partial t} - \nabla \cdot \boldsymbol{\sigma} = -(\rho_e \mathbf{E} + \mathbf{J} \times \mathbf{B})$$

## 3. Energy Conditions
- **Weak Energy Condition (WEC)**: $T_{\mu\nu} v^\mu v^\nu \ge 0$ for all timelike vectors $v^\mu$.
- **Null Energy Condition (NEC)**: $T_{\mu\nu} k^\mu k^\nu \ge 0$ for all null vectors $k^\mu$.
- **Dominant Energy Condition (DEC)**: $-T^\mu{}_\nu v^\nu$ is future-directed causal vector for timelike $v^\mu$, ensuring energy flux does not travel faster than light ($|\mathbf{S}|/c \le u$).
- **Strong Energy Condition (SEC)**: $\left( T_{\mu\nu} - \frac{1}{2} g_{\mu\nu} T^\lambda{}_\lambda \right) v^\mu v^\nu \ge 0$.
