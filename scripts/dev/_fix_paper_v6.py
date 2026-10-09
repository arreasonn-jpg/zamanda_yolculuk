"""Add 6 figures to main.tex at logical positions."""
from pathlib import Path

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# Guard: if figures already added, skip
if "\\includegraphics" in c:
    print("[SKIP] Figures already present")
    exit(0)

# Figure blocks to insert
fig1 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig1_unified_mechanisms.png}
\caption{Feasibility gap across six mechanisms. Bars show
$\log_{10}(\text{required}/\text{device})$; green = physically
allowed, red = excluded by Hawking (1992).}
\label{fig:unified}
\end{figure}
"""

fig2 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig2_em_scaling.png}
\caption{EM-induced metric perturbation $h_{00}$ vs.\ drive current $I$.
Reference lines mark lightning, lab-pulsed magnets, astrophysical
plasmas, and magnetar-scale currents. The $1$~s/yr target requires
$I \sim 10^{22}$~A.}
\label{fig:scaling}
\end{figure}
"""

fig3 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig3_best_case_tiers.png}
\caption{Best-case engineering tiers (baseline copper, cooled copper,
superconducting HTS with multi-stage array). Even the most aggressive
scenario remains $10^{30}$ below the target.}
\label{fig:tiers}
\end{figure}
"""

fig4 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig4_monte_carlo.png}
\caption{Monte Carlo propagation of parameter uncertainties
($N=50{,}000$). Relative standard deviation on $h_{00}$ is 14.6\%.}
\label{fig:mc}
\end{figure}
"""

fig5 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig5_sensitivity_tornado.png}
\caption{Sensitivity indices $S_i=(dh/h)/(dp/p)$ for device parameters.
Current and coil radius dominate.}
\label{fig:tornado}
\end{figure}
"""

fig6 = r"""
\begin{figure}[t]
\centering
\includegraphics[width=\columnwidth]{figures/fig6_backpack_layout.png}
\caption{Backpack component layout within the $45\times 32\times 15$~cm
envelope. Point size scales with component mass.}
\label{fig:backpack}
\end{figure}
"""

# Insert Fig 1 (unified) at start of Discussion
marker = "\\section{Discussion}"
c = c.replace(marker, fig1 + "\n" + marker, 1)

# Insert Fig 2 (scaling) after the paragraph that references fig:scaling
marker2 = "Fig.~\\ref{fig:scaling} shows"
if marker2 in c:
    # Find end of paragraph (next blank line after marker2)
    idx = c.find(marker2)
    end_para = c.find("\n\n", idx)
    if end_para > 0:
        c = c[:end_para] + "\n" + fig2 + c[end_para:]

# Insert Fig 3 (tiers) after best-case paragraph
marker3 = "a $Q$ factor of $10^{34}$."
if marker3 in c:
    idx = c.find(marker3)
    end_para = c.find("\n\n", idx)
    if end_para > 0:
        c = c[:end_para] + "\n" + fig3 + c[end_para:]

# Insert Fig 4 (MC) + Fig 5 (tornado) at end of Section 3.4
marker4 = "for the region diameter."
if marker4 in c:
    idx = c.find(marker4)
    end_para = c.find("\n\n", idx)
    if end_para > 0:
        c = c[:end_para] + "\n" + fig4 + "\n" + fig5 + c[end_para:]

# Insert Fig 6 (backpack) in Method section
marker5 = "\\subsection{Validation}"
if marker5 in c:
    c = c.replace(marker5, fig6 + "\n" + marker5, 1)

P.write_text(c, encoding="utf-8")

# Report
import re
figs = len(re.findall(r"\\includegraphics", c))
print(f"[OK] Inserted {figs} figure(s)")
print()
labels = re.findall(r"\\label\{(fig:[^}]+)\}", c)
print("Figure labels:", labels)