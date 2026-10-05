"""Fix main.tex: email, duplicate sections, Table 2 overflow."""
from pathlib import Path
import re

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# ===== DIAGNOSTIC =====
print("=== BEFORE ===")
sections = re.findall(r"\\section\*?\{([^}]+)\}", c)
for i, s in enumerate(sections, 1):
    print(f"  {i}. {s}")
disc_count = len(re.findall(r"\\section\{Discussion\}", c))
print(f"\\section{{Discussion}} count: {disc_count}")
print()

# ===== FIX 1: Email =====
c = c.replace("a.rreasonn@gmail.com", "erdem.esa.71@gmail.com")

# ===== FIX 2: Duplicate Discussion/Conclusion =====
disc_positions = [m.start() for m in re.finditer(r"\\section\{Discussion\}", c)]
if len(disc_positions) > 1:
    print(f"Removing duplicate Discussion at position {disc_positions[1]}...")
    bibstyle_pos = c.find("\\bibliographystyle")
    if bibstyle_pos < 0:
        bibstyle_pos = len(c)
    c = c[:disc_positions[1]] + c[bibstyle_pos:]

# ===== FIX 3: Table 2 (unified feasibility) - make it fit =====
# Find the FIRST table* environment (Table 2)
start = c.find("\\begin{table*}")
end = c.find("\\end{table*}")
if start > 0 and end > 0:
    new_table = r"""\begin{table*}[t]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Unified feasibility across six mechanisms.
Gap $=$ required resource $/$ device-scale resource.}
\label{tab:unified}
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
\end{table*}
"""
    c = c[:start] + new_table + c[end + len("\\end{table*}"):]
    print("[OK] Table 2 shortened to 5 columns + footnotesize")

# ===== FIX 4: Section reference in Intro =====
c = c.replace("Section 7 summarizes our conclusions.",
              "Section \\ref{sec:conclusion} summarizes our conclusions.")

# ===== SAVE =====
P.write_text(c, encoding="utf-8")

# ===== FINAL DIAGNOSTIC =====
print()
print("=== AFTER ===")
sections = re.findall(r"\\section\*?\{([^}]+)\}", c)
for i, s in enumerate(sections, 1):
    print(f"  {i}. {s}")
disc_count2 = len(re.findall(r"\\section\{Discussion\}", c))
print(f"\\section{{Discussion}} count: {disc_count2}")
print()
print("[OK] main.tex updated")