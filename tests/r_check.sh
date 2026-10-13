#!/usr/bin/env bash
# Run R/downstream_analysis.R on synthetic data with 6 known enriched genera and check
# that the DA consensus recovers exactly those genera.
# Usage: bash tests/r_check.sh WORKDIR     (WORKDIR must not exist yet; needs Rscript + packages)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/r_check.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK/export"

Rscript "${REPO}/tests/make_synthetic.R" "$WORK/export"
Rscript "${REPO}/R/downstream_analysis.R" "$WORK/export" "$WORK/out" treatment > "$WORK/run.log" 2>&1 \
    || { tail -20 "$WORK/run.log"; exit 1; }

found="$(awk -F, 'NR > 1 { gsub(/"/, "", $1); sub(/.*g__/, "", $1); print $1 }' "$WORK/out/da_consensus.csv" | sort | tr '\n' ' ')"
expected="Genus01 Genus02 Genus03 Genus04 Genus05 Genus06 "
echo "consensus: ${found}"
if [[ "$found" == "$expected" ]]; then echo "PASS"; else echo "FAIL (expected: ${expected})"; exit 1; fi
