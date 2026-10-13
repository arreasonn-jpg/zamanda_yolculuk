"""v7.4-alpha2: 4-regime diagram for the paper."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = np.load("results/validated/v7_4_sub_qft.npz")
A = data["A"]; rho_min = data["rho_min"]; ratio = data["ratio"]

fig, ax = plt.subplots(1, 1, figsize=(9, 5))

# Plot points colored by regime
for a, r, ra in zip(A, rho_min, ratio):
    if r == 0.0:
        color, label = "gray", "Deactivated/No NEC" if a == 0.5 else None
    elif ra < 1.0:
        color, label = "green", "NEC + QFT-OK" if a == 1.4 else None
    else:
        color, label = "red", "NEC + QFT-violating" if a == 2.0 else None
    ax.scatter(a, r, color=color, s=80, zorder=3,
                label=label, edgecolors='k', linewidths=0.5)

# Vertical threshold lines
ax.axvline(1.2, color='gray', ls=':', lw=1, alpha=0.7)
ax.axvline(1.35, color='orange', ls=':', lw=1, alpha=0.7)
ax.axvline(1.87, color='red', ls='--', lw=1.5)

# Annotations
ax.text(1.2, 0.3, r'$A_{\rm act} \approx 1.2$', rotation=90,
        fontsize=9, color='gray', va='bottom')
ax.text(1.35, 0.3, r'$A_{\rm NEC} \approx 1.35$', rotation=90,
        fontsize=9, color='orange', va='bottom')
ax.text(1.87, 0.3, r'$A_{\rm QFT} \approx 1.87$', rotation=90,
        fontsize=9, color='red', va='bottom')

ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel(r'$A_{\rm fuel}$', fontsize=11)
ax.set_ylabel(r'$\rho_{\min}$', fontsize=11)
ax.set_title('PTUM 4-regime structure', fontsize=12)
ax.grid(True, alpha=0.3)
ax.legend(loc='lower left', fontsize=9)

# Regime labels at top
ax.text(0.85, ax.get_ylim()[0]*0.15, 'Deactivated',
        ha='center', fontsize=9, style='italic')
ax.text(1.55, ax.get_ylim()[0]*0.15, 'NEC + QFT-OK',
        ha='center', fontsize=9, style='italic', color='green')
ax.text(2.6, ax.get_ylim()[0]*0.15, 'NEC + QFT-violating',
        ha='center', fontsize=9, style='italic', color='red')

plt.tight_layout()
plt.savefig("paper/figures/fig9_regimes.png", dpi=110)
print("Saved: paper/figures/fig9_regimes.png")
