"""Fix main.tex: email + duplicate sections + Table 2 width."""
from pathlib import Path
import re

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# ===== DIAGNOSTIC BEFORE =====
print("=== BEFORE ===")
for i, m in enumerate(re.finditer(r"\\section\*?\{([^}]+)\}", c), 1):
    print(f"  {i}. {m.group(1)}")

# ===== FIX 1: Email =====
c = c.replace("a.rreasonn@gmail.com", "erdem.esa.71@gmail.com")

# ===== FIX 2: Nuke all Discussion/Conclusion/Ack/Data sections =====
# Find FIRST \section{Discussion}
start = c.find("\\section{Discussion}")
if start < 0:
    print("ERROR: no Discussion section")
    exit(1)

# Find \bibliographystyle
end = c.find("\\bibliographystyle")
if end < 0:
    print("ERROR: no bibliographystyle")
    exit(1)

# ===== CLEAN REPLACEMENT =====
clean = r"""\section{Discussion}
\label{sec:discussion}

\subsection{Physical interpretation}
\label{subsec:interp}

The systematic failure of all six mechanisms to reach device-scale
feasibility is not accidental. Each mechanism requires either:
(a) an energy density at or above the QCD scale (EM, wormholes,
Alcubierre), (b) an angular momentum comparable to astrophysical
bodies (Kerr, Tipler, G\"odel), or (c) infinite spatial extent
(Tipler). The EM case is particularly instructive: even though
electromagnetic fields are the most experimentally accessible
spacetime-curving source available to a laboratory, their energy
density at achievable field strengths is too low by $\sim 40$
orders of magnitude.

\subsection{Comparison with prior work}
\label{subsec:prior}

Our EM bound is consistent with order-of-magnitude estimates
scattered through the literature, but the present work provides,
to our knowledge, the first full first-principles pipeline for
a specific device geometry. The Tipler and Morris--Thorne bounds
match the original estimates in \citet{Tipler1974} and
\citet{MorrisThorne1988}. The Alcubierre estimate is somewhat
larger than the popular ``Jupiter mass'' figure often quoted in
secondary sources \citep{Alcubierre1994}, because we evaluate the
energy at $R=100$~m and $v=c$ rather than the smaller radii used
in some early optimizations.

\subsection{Implications}
\label{subsec:impl}

Three implications follow. First, no electromagnetic configuration
realizable in a backpack-sized device can produce a macroscopic
temporal effect. Second, the same conclusion holds for every
classical GR mechanism that could in principle generate a CTC.
Third, if a device-scale time-dilation experiment is ever to be
attempted, it must invoke physics beyond classical GR.

\subsection{Limitations}
\label{subsec:limits}

Our analysis is confined to classical GR coupled to classical EM.
We do not consider: (i) semiclassical effects near Planck-scale
field strengths; (ii) non-minimal coupling between EM and gravity;
(iii) modified-gravity scenarios; (iv) quantum-gravity corrections.
Any of these could in principle shift the feasibility threshold,
though none is currently supported by observation. We also assume
ideal lossless coils; real coils will show slightly larger $h_{00}$
per unit current, but not by a factor that closes the $10^{20}$ gap.

\subsection{Future work}
\label{subsec:future}

Extensions of this work could include: a full 3+1D general-
relativistic simulation of the strong-field regime; inclusion of
semiclassical corrections; systematic optimization of coil
geometry; and a public release of the pipeline for independent
verification.

\section{Conclusion}
\label{sec:conclusion}

We have presented a unified feasibility analysis of six
physically-motivated mechanisms for producing macroscopic time
dilation at device scale. The electromagnetic case is treated
with a first-principles numerical pipeline validated to
sub-percent accuracy; the other five are bounded by published
closed-form expressions. All six mechanisms require resources
between $10^{20}$ and $10^{49}$ times beyond device scale.
The closest physically-allowed mechanism (electromagnetic) still
requires $\sim 10^{22}$~A, exceeding magnetar-scale currents by
a factor that cannot be closed by any current or foreseeable
engineering improvement. We conclude that no
classical-GR-consistent mechanism for macroscopic time dilation
is accessible at device scale.

\section*{Acknowledgments}
The author thanks the open-source scientific Python community for
the tools (\texttt{numpy}, \texttt{scipy}, \texttt{matplotlib})
that made this work possible.

\section*{Data availability}
All code, figures, and validation data are available at the
project repository.

"""

# Replace from Discussion start to bibliographystyle
c = c[:start] + clean + c[end:]

# ===== FIX 3: Section ref in Intro =====
c = c.replace("Section 7 summarizes our conclusions",
              "Section~\\ref{sec:conclusion} summarizes our conclusions")

# ===== FIX 4: Table 2 — already footnotesize; add resizebox if missing =====
if "resizebox" not in c:
    # Find first table* (unified feasibility)
    ts = c.find("\\begin{table*}")
    te = c.find("\\end{table*}")
    if ts > 0 and te > 0:
        c = (c[:ts]
             + "\\begin{table*}[t]\n\\centering\n\\small\n\\resizebox{\\textwidth}{!}{%\n"
             + c[ts + len("\\begin{table*}[t]"):te]
             + "}\n\\end{table*}"
             + c[te + len("\\end{table*}"):])

P.write_text(c, encoding="utf-8")

# ===== DIAGNOSTIC AFTER =====
print()
print("=== AFTER ===")
for i, m in enumerate(re.finditer(r"\\section\*?\{([^}]+)\}", c), 1):
    print(f"  {i}. {m.group(1)}")
print()
print("[OK] main.tex rewritten")