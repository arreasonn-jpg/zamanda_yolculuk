"""Small domain so wave hits BC within 50 steps."""
import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta3 import rk4_step_beta3
from zt005.bssn_rk4_beta2 import rk4_step_beta2

N = 9; L = 1.0; h = L/(N-1); dt = 1e-3
def make():
    s = minkowski_state(N)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = 1e-3 * np.exp(-(X**2 + Y**2 + Z**2) / 0.05)
    return s

s_no = make()
s_bc = make()
for step in range(100):
    s_no = rk4_step_beta2(s_no, h=h, dt=dt)
    s_bc = rk4_step_beta3(s_bc, h=h, dt=dt)

phi_no = np.max(np.abs(s_no.phi))
phi_bc = np.max(np.abs(s_bc.phi))
print(f"without BC: max|phi| = {phi_no:.4e}")
print(f"with BC:    max|phi| = {phi_bc:.4e}")
print(f"ratio: {phi_bc/phi_no:.4f}")
print(f"BC effect: {(1-phi_bc/phi_no)*100:.2f}%")
