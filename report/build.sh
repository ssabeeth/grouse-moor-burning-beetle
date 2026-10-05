#!/bin/bash
# Build both reports with xelatex: report_v6.pdf (heather management, presented as "Version 2")
# and beetle_report.pdf (heather beetle damage). Both share style.tex.
# Run from claude_code_kit/ first:
#   python report/make_figures.py && python report/make_beetle_figures.py
set -e
cd "$(dirname "$0")"
for doc in report_v6 beetle_report; do
  latexmk -xelatex -interaction=nonstopmode -halt-on-error -quiet "$doc.tex"
  latexmk -c -quiet "$doc.tex" && rm -f "$doc.xdv"
done
