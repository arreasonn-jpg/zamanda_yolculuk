"""Fix main.tex: resize Table 1 (tab:em_results)."""
from pathlib import Path

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# Find the FIRST table environment that contains \label{tab:em_results}
# (Table 2 uses table* so this one is table)
label_pos = c.find("\\label{tab:em_results}")
if label_pos < 0:
    print("ERROR: tab:em_results label not found")
    exit(1)

# Find enclosing \begin{table}...\end{table}
start = c.rfind("\\begin{table}", 0, label_pos)
end = c.find("\\end{table}", label_pos)

if start < 0 or end < 0:
    print("ERROR: table environment not found")
    exit(1)

clean_table = r"""\begin{table}[t]
\centering
\small
\caption{EM-induced metric perturbation at device center
(baseline 100~A, 1~MHz).}
\label{tab:em_results}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{lcc}
\toprule
Geometry & $h_{00}$ (center) & $|h_{0i}|$ (center) \\
\midrule
G1 (solenoid)          & $5.94\times 10^{-48}$ & $<10^{-50}$ \\
G2 (toroid)            & $1.24\times 10^{-49}$ & $<10^{-50}$ \\
G3 (counter-rotating)  & $2.03\times 10^{-49}$ & $<10^{-51}$ \\
G4 (hybrid)            & $9.36\times 10^{-49}$ & $<10^{-50}$ \\
\bottomrule
\end{tabular}
}
\end{table}
"""

c = c[:start] + clean_table + c[end + len("\\end{table}"):]
P.write_text(c, encoding="utf-8")
print("[OK] Table 1 rewritten with resizebox")