"""Fix: add lmodern + disable microtype expansion."""
from pathlib import Path

P = Path("paper/main.tex")
c = P.read_text(encoding="utf-8")

# Add lmodern after fontenc
if "lmodern" not in c:
    c = c.replace(
        "\\usepackage[T1]{fontenc}",
        "\\usepackage[T1]{fontenc}\n\\usepackage{lmodern}"
    )
    print("[OK] Added lmodern")

# Replace microtype with expansion disabled
c = c.replace(
    "\\usepackage{microtype}",
    "\\usepackage[protrusion=true,expansion=false]{microtype}"
)
print("[OK] microtype: expansion disabled")

P.write_text(c, encoding="utf-8")
print()
print("Now rebuild FROM SCRATCH (delete .aux first):")
print("  cd paper")
print("  Remove-Item main.aux,main.bbl,main.blg,main.out,main.log -ErrorAction SilentlyContinue")
print("  pdflatex -interaction=nonstopmode main.tex")
print("  bibtex main")
print("  pdflatex -interaction=nonstopmode main.tex")
print("  pdflatex -interaction=nonstopmode main.tex")