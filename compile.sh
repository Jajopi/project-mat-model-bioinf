#!/usr/bin/env bash
set -euo pipefail
input="${1:-results.md}"
output="${2:-${input%.md}.pdf}"
pandoc "$input" --template=template.latex --pdf-engine=pdflatex -o "$output"
