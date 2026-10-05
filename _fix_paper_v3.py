"""Fix main.tex: rewrite Table 2 with proper label placement."""
from pathlib import Path

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# Find table* block (unified feasibility = Table 2)
ts = c.find("\\begin{table*}")
te = c.find("\\end{table*}")

if ts < 0 or te < 0:
    print("ERROR: table* block not found")
    exit(1)

clean_table = r"""\begin{table*}[t]
\centering
\caption{Unified feasibility across six mechanisms.
Gap $=$ required resource $/$ device-scale value.}
\label{tab:unified}
\resizebox{0.95\textwidth}{!}{%
\begin{tabular}{llrrr}
\toprule
\# & Mechanism & Parameter & Required & Gap \\
\midrule
1 & EM field          & $I$ (A)                 & $10^{22}$ & $10^{20}$ \\
2 & Kerr              & $J$ (kg\,m$^2$/s)       & $10^{27}$ & $10^{48}$ \\
3 & Tipler            & $\rho$ (kg/m$^3$)       & $10^{23}$ & $10^{20}$ \\
4 & Morris--Thorne    & $|\rho|$ (J/m$^3$)      & $10^{49}$ & $10^{41}$ \\
5 & Alcubierre        & $|\rho|$ (J/m$^3$)      & $10^{57}$ & $10^{49}$ \\
6 & G\"odel           & $\rho$ (kg/m$^3$)       & $10^{26}$ & $10^{23}$ \\
\bottomrule
\end{tabular}
}
\end{table*}
"""

c = c[:ts] + clean_table + c[te + len("\\end{table*}"):]
P.write_text(c, encoding="utf-8")

# Diagnostics
import re
labels = re.findall(r"\\label\{([^}]+)\}", c)
refs = re.findall(r"\\ref\{([^}]+)\}", c)
print("=== LABELS ===")
for l in labels:
    print(f"  {l}")
print()
print("=== REFS ===")
for r in refs:
    print(f"  {r}")
print()
print("=== ORPHAN REFS (ref without label) ===")
for r in refs:
    if r not in labels:
        print(f"  [X] {r}")
print()
print("[OK] main.tex Table 2 rewritten")