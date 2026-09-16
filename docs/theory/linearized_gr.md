# Linearized General Relativity and Retarded Green Potentials

## 1. Linearized Einstein Equations
In weak-field gravity where $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ with $|h_{\mu\nu}| \ll 1$, defining the trace-reversed perturbation:
$$\bar{h}_{\mu\nu} = h_{\mu\nu} - \frac{1}{2} \eta_{\mu\nu} h^\lambda{}_\lambda$$
and adopting the Lorenz / harmonic gauge condition:
$$\partial_\mu \bar{h}^{\mu\nu} = 0$$
reduces the Einstein equations to the wave equation:
$$\Box \bar{h}_{\mu\nu} = -\frac{16\pi G}{c^4} T_{\mu\nu}$$

## 2. Retarded Potential Solution
The physical causal solution satisfying outgoing radiation boundary conditions is:
$$\bar{h}_{\mu\nu}(t, \mathbf{x}) = \frac{4G}{c^4} \int \frac{T_{\mu\nu}(t - |\mathbf{x} - \mathbf{x}'|/c, \mathbf{x}')}{|\mathbf{x} - \mathbf{x}'|} d^3x'$$

In the static / quasistatic limit ($f \to 0$ or retardation $|\mathbf{x}-\mathbf{x}'|/c \ll 1/f$):
$$\bar{h}_{\mu\nu}(\mathbf{x}) = \frac{4G}{c^4} \int \frac{T_{\mu\nu}(\mathbf{x}')}{|\mathbf{x} - \mathbf{x}'|} d^3x'$$
Physical metric perturbation $h_{\mu\nu}$ is reconstructed via:
$$h_{\mu\nu} = \bar{h}_{\mu\nu} - \frac{1}{2} \eta_{\mu\nu} \bar{h}^\lambda{}_\lambda$$
