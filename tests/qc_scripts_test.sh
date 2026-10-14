#!/usr/bin/env bash
# Unit tests for pipeline/denoise_check.py and pipeline/choose_depth.py on small inputs with
# hand-computed answers. Needs python3 only.
# Usage: bash tests/qc_scripts_test.sh WORKDIR     (WORKDIR must not exist yet)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/qc_scripts_test.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
status=0
same() {  # same DESCRIPTION GOT WANT
    if [[ "$2" == "$3" ]]; then echo "PASS  $1"; else echo "FAIL  $1: got '$2', want '$3'"; status=1; fi
}

# ---- denoise_check.py: 10 typical samples, one poor filter, one poor merge, one shallow
{
    printf 'sample-id\tinput\tfiltered\tpercentage of input passed filter\tdenoised\tmerged\tpercentage of input merged\tnon-chimeric\tpercentage of input non-chimeric\n'
    printf '#q2:types\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\n'
    for i in 1 2 3 4 5 6 7 8 9 10; do
        printf 'T%s\t10000\t%s\t0\t%s\t%s\t0\t%s\t0\n' "$i" $((7800 + i * 10)) $((7600 + i * 10)) $((7350 + i * 10)) $((7340 + i * 10))
    done
    printf 'BADFILT\t10000\t5000\t0\t4900\t4750\t0\t4740\t0\n'     # 50% passed
    printf 'BADMERGE\t10000\t7800\t0\t7600\t5700\t0\t5690\t0\n'    # 75% of denoised merged
    printf 'SHALLOW\t1000\t780\t0\t760\t735\t0\t734\t0\n'          # typical rates, 734 reads
} > "$WORK/stats.tsv"
python3 "$REPO/pipeline/denoise_check.py" --stats "$WORK/stats.tsv" --out "$WORK/flags.tsv" \
    --min-depth 1000 > "$WORK/check.log"
flag() { awk -F'\t' -v s="$1" '$1 == s { print $7 }' "$WORK/flags.tsv"; }
same "poor filter pass flagged"              "$(flag BADFILT)" low_filter_pass
same "poor merging flagged"                  "$(flag BADMERGE)" low_merge
same "shallow sample flagged only for depth" "$(flag SHALLOW)" low_depth
same "typical samples not flagged"           "$(flag T1)$(flag T5)$(flag T10)" ""
same "no chimera flags when chimeras are ~0" "$(grep -c high_chimera "$WORK/flags.tsv" || true)" 0

# ---- choose_depth.py: two groups; main = smallest sample; sensitivity keeps >= 80% per group
printf 'sample\treads\nA1\t1000\nA2\t1500\nA3\t2000\nA4\t2500\nA5\t3000\nB1\t1100\nB2\t5000\nB3\t6000\nB4\t7000\nB5\t8000\n' > "$WORK/depths.tsv"
printf 'sample-id\tsite\tlayer\n#q2:types\tcategorical\tcategorical\nA1\tA\ttop\nA2\tA\ttop\nA3\tA\ttop\nA4\tA\ttop\nA5\tA\ttop\nB1\tB\ttop\nB2\tB\ttop\nB3\tB\ttop\nB4\tB\ttop\nB5\tB\ttop\n' > "$WORK/meta.tsv"
got="$(python3 "$REPO/pipeline/choose_depth.py" --depths "$WORK/depths.tsv" --metadata "$WORK/meta.tsv" \
    --groups "site layer" --out "$WORK/dc.tsv" --groups-out "$WORK/dg.tsv")"
# group A sorted desc 3000 2500 2000 1500 1000: 4 of 5 must stay -> 1500; group B -> 5000; min = 1500
same "main 1000, sensitivity 1500 (limited by site A)" "$got" "1000 1500"
same "limiting group recorded" "$(awk -F'\t' '$1 == "limiting_groups" { print $2 }' "$WORK/dc.tsv")" "A / top"
same "dropped samples listed per group" "$(awk -F'\t' '$1 == "A / top" { print $2, $3, $4, $5 }' "$WORK/dg.tsv")" "5 5 4 A1"
same "no groups: main depth only" \
     "$(python3 "$REPO/pipeline/choose_depth.py" --depths "$WORK/depths.tsv" --out "$WORK/dc2.tsv")" 1000
printf 'sample-id\tsite\nA1\tA\n' > "$WORK/meta_short.tsv"
rc=0
python3 "$REPO/pipeline/choose_depth.py" --depths "$WORK/depths.tsv" --metadata "$WORK/meta_short.tsv" \
    --groups site --out "$WORK/x.tsv" 2> "$WORK/err.txt" || rc=$?
same "samples missing from metadata stop with a message" \
     "$rc $(grep -c "missing from the metadata" "$WORK/err.txt")" "1 1"

exit "$status"
