#!/usr/bin/env bash
set -euo pipefail
input="${1:-README.md}"
output="${2:-results.pdf}"
pandoc "$input" --template=template.latex --pdf-engine=pdflatex -o "$output"
