"""Build paper/archive/arxiv_submission_v4.zip from the live paper/arxiv_v3
sources (snapshot dated 2026-10-20, repo state v13.0.0).

Contents mirror the v3 zip convention: arxiv_v3/ + tex + refs.bib + figures/.
Regenerate after further updates (v13.1.0 will regenerate with the 4x
subcycle result written into the bssn sections).
"""
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "paper" / "arxiv_v3"
OUT = ROOT / "paper" / "archive" / "arxiv_submission_v4.zip"

KLUDGE_SAFE = ("main_en.tex", "main_tr.tex", "bssn_section_en.tex",
               "bssn_section_tr.tex", "refs.bib")

files = [f for f in KLUDGE_SAFE if (SRC / f).exists()]
figs = sorted((SRC / "figures").glob("*.png"))
assert figs, "no figures found"
assert len(files) == len(KLUDGE_SAFE), f"missing tex files: {files}"

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for f in files:
        z.write(SRC / f, f"arxiv_v3/{f}")
    for f in figs:
        z.write(f, f"arxiv_v3/figures/{f.name}")

print(f"wrote {OUT}  ({len(files)} tex + {len(figs)} figures)")
