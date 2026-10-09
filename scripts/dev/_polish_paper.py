"""Polish main.tex: microtype + authblk + Turkish abstract + appendix."""
from pathlib import Path

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# Guard: skip if already polished
if "microtype" in c:
    print("[SKIP] Already polished")
    exit(0)

# ===== 1. Add packages =====
c = c.replace(
    "\\usepackage{natbib}",
    "\\usepackage{natbib}\n"
    "\\usepackage{microtype}\n"
    "\\usepackage{babel}\n"
    "\\usepackage{authblk}"
)

# ===== 2. Replace author block with authblk format =====
old_author = r"""\author{
  Erdem Esa\thanks{Corresponding author.} \\
  \small Independent Researcher \\
  \small \texttt{erdem.esa.71@gmail.com}
}"""

new_author = r"""\author[1]{Erdem Esa\thanks{Corresponding author:
\texttt{erdem.esa.71@gmail.com}}}
\affil[1]{Independent Researcher}"""

if old_author in c:
    c = c.replace(old_author, new_author)
    print("[OK] Author block -> authblk")
else:
    print("[WARN] Author block not matched, skipping")

# ===== 3. Add Turkish abstract after English abstract =====
tr_abstract_marker = "\\end{abstract}"
tr_abstract = r"""
\begin{otherlanguage}{turkish}
\renewcommand{\abstractname}{T\"urk\c{c}e \"Ozet}
\begin{abstract}
Bu \c{c}al\i{}\c{s}mada, laboratuvar \"ol\c{c}e\u{g}inde makroskopik
zaman geni\c{s}lemesi veya kapal\i{} zaman benzeri e\u{g}riler
(CTC) \"uretmek i\c{c}in \"onerilen alt\i{} fiziksel mekanizma
birle\c{s}ik bir say\i{}sal ve analitik \c{c}er\c{c}evede
incelenmi\c{s}tir: (i) elektromanyetik alan konfig\"urasyonlar\i{},
(ii) Kerr \c{c}er\c{c}eve s\"ur\"uklenmesi, (iii) Tipler d\"onen
silindirleri, (iv) Morris--Thorne solucan delikleri, (v) Alcubierre
b\"uk\"um s\"ur\"uc\"uleri ve (vi) G\"odel d\"onen evreni.
Elektromanyetik durum i\c{c}in, sonlu-telli Biot--Savart kaynak
modeli, tam elektromanyetik enerji-momentum tens\"or\"u ve
lineerle\c{s}tirilmi\c{s} Einstein denklemleri i\c{c}in eksenel
simetrik sonlu fark Poisson \c{c}\"oz\"uc\"us\"un\"u birle\c{s}tiren
ilk-ilkelerden t\"uretilmi\c{s} bir boru hatt\i{} geli\c{s}tirdik.
Di\u{g}er be\c{s} mekanizma i\c{c}in kapal\i{} formda
fizibilite s\i{}n\i{}rlar\i{} t\"urettik. Alt\i{} mekanizman\i{}n
tamam\i{}, cihaz \"ol\c{c}e\u{g}inin $10^{20}$ ile $10^{49}$ kat\i{}
aras\i{}nda kaynak gerektirmektedir. Fiziksel olarak izin verilen
en yak\i{}n mekanizma (elektromanyetik), h\^al\^a magnetar
\"ol\c{c}e\u{g}indeki ak\i{}mlar\i{} a\c{s}an $\sim 10^{22}$~A
ak\i{}m gerektirmektedir. Hi\c{c}bir m\"uhendislik iyile\c{s}tirmesi
bu bo\c{s}lu\u{g}u kapatamaz. Sonu\c{c} olarak, klasik genel
g\"orelilik ile tutarl\i{} hi\c{c}bir mekanizman\i{}n cihaz
\"ol\c{c}e\u{g}inde makroskopik zaman geni\c{s}lemesi i\c{c}in
eri\c{s}ilebilir olmad\i{}\u{g}\i{}n\i{} g\"osteriyoruz.
\end{abstract}
\end{otherlanguage}"""

if tr_abstract_marker in c and "otherlanguage" not in c:
    c = c.replace(tr_abstract_marker, tr_abstract_marker + tr_abstract, 1)
    print("[OK] Turkish abstract added")
else:
    print("[WARN] Turkish abstract not added (marker missing or already present)")

# ===== 4. Add appendix before \end{document} =====
appendix_block = r"""
\appendix
\section{Reproducibility}
\label{app:repro}

All numerical results in this paper were produced by the open-source
Python pipeline described in Section~\ref{sec:method}. The complete
source code, validation scripts, and generated figures are available
at the project repository.

\paragraph{Software dependencies.}
\texttt{numpy}~$\geq$~2.0, \texttt{scipy}~$\geq$~1.10,
\texttt{matplotlib}~$\geq$~3.8, \texttt{pytest}~$\geq$~8.

\paragraph{Reproduction.}
Run \texttt{python -m pytest} for the unit-test suite (60+ tests)
and \texttt{python -m zt005.run\_zt006\_3} for the EM metric
pipeline. Numerical results in Tables~\ref{tab:em_results}
and~\ref{tab:unified} are reproduced by
\texttt{run\_zt007\_b7.py}.

\paragraph{Grid convergence.}
The grid convergence study reported in Section~\ref{subsec:validation}
uses the script \texttt{run\_zt006\_2\_23.py} with \texttt{--mode full}.

"""

end_doc = "\\end{document}"
if "\\appendix" not in c:
    c = c.replace(end_doc, appendix_block + end_doc, 1)
    print("[OK] Appendix added")

P.write_text(c, encoding="utf-8")
print()
print("[OK] main.tex polished")
print()
print("Recompile:")
print("  cd paper")
print("  pdflatex -interaction=nonstopmode main.tex")
print("  bibtex main")
print("  pdflatex -interaction=nonstopmode main.tex")
print("  pdflatex -interaction=nonstopmode main.tex")