#!/usr/bin/env bash
# Run R/design_analysis.R on a synthetic split-plot experiment with known effects
# (tests/make_synthetic_splitplot.R: 18 plots x 3 exclosure subplots x 4 depths) and check
# the strata, degrees of freedom, effects and DA consensus. Needs Rscript with vegan,
# lme4, lmerTest, ANCOMBC, maaslin3, ALDEx2, ggplot2; and python3 (depth-rule parity).
# Takes about 10 minutes.
# Usage: bash tests/design_check.sh WORKDIR     (WORKDIR must not exist yet)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/design_check.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
status=0
same() {  # same DESCRIPTION GOT WANT
    if [[ "$2" == "$3" ]]; then echo "PASS  $1"; else echo "FAIL  $1: got '$2', want '$3'"; status=1; fi
}
within() {  # within DESCRIPTION VALUE LOW HIGH
    if awk -v v="$2" -v lo="$3" -v hi="$4" 'BEGIN { exit !(v != "" && v >= lo && v <= hi) }'; then
        echo "PASS  $1 ($2)"
    else echo "FAIL  $1: $2 not in [$3, $4]"; status=1; fi
}
# col FILE COLUMN [FILTER_COLUMN VALUE ...] -> values of COLUMN in rows matching all filters
col() {
    local f="$1" c="$2"; shift 2
    awk -F'\t' -v c="$c" -v filt="$*" '
        NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; n = split(filt, fl, " "); next }
        { ok = 1; for (k = 1; k < n; k += 2) if ($(h[fl[k]]) != fl[k + 1]) ok = 0
          if (ok) print $(h[c]) }' "$f"
}
genera() { col "$1" feature coefficient "$2" | sed 's/.*g__//' | sort | tr '\n' ' ' | sed 's/ $//'; }

Rscript "$REPO/tests/make_synthetic_splitplot.R" "$WORK/data" > /dev/null
design() {  # design OUT_FILE TEMPLATE [sed expressions...]
    local outf="$1" tmpl="$2"; shift 2
    sed -e "s#^export_dir\t.*#export_dir\t$WORK/data#" -e 's/^permutations\t.*/permutations\t199/' \
        -e 's/^rarefy_iter\t.*/rarefy_iter\t10/' -e 's/^cores\t.*/cores\t2/' "$@" \
        "$REPO/R/designs/$tmpl" > "$outf"
}

# ---- A. Manuscript-2 design: all depths, split-plot with depth as a third stratum ----
design "$WORK/ms2.tsv" ms2_community_all_depths.tsv
Rscript "$REPO/R/design_analysis.R" "$WORK/ms2.tsv" "$WORK/ms2" > "$WORK/ms2.log" 2>&1 \
    || { echo "FAIL  design_analysis.R (ms2) failed: $(tail -5 "$WORK/ms2.log")"; exit 1; }
A="$WORK/ms2"
same "strata: nutrient/round -> plot, exclosure -> subplot, depth -> sample" \
     "$(for t in nutrient round exclosure depth_layer nutrient:exclosure:depth_layer; do col "$A/design_strata.tsv" stratum term "$t"; done | tr '\n' ' ')" \
     "plot plot subplot sample sample "
same "samples below MIN_SAMPLE_DEPTH removed" "$(col "$A/samples_removed.tsv" sample | sort | tr '\n' ' ')" \
     "GA1UC7S04 GA2FB4S08 "
same "random intercepts for plot and subplot" "$(col "$A/design_used.tsv" value setting random_effects)" \
     "(1 | plot) + (1 | subplot)"
# split-plot denominator df: plots 18 - 6 cells = 12; subplots 54 - 18 - 12 = 24
within "alpha LMM: nutrient tested on plot df" "$(col "$A/alpha_lmm.tsv" DenDF metric Shannon term nutrient)" 11 13
within "alpha LMM: exclosure tested on subplot df" "$(col "$A/alpha_lmm.tsv" DenDF metric Shannon term exclosure)" 22 26
within "alpha LMM: depth tested on sample df" "$(col "$A/alpha_lmm.tsv" DenDF metric Shannon term depth_layer)" 90 130
P="$A/permanova_hierarchical.tsv"
same "PERMANOVA reports all 15 terms per distance" "$(col "$P" term distance bray_curtis | wc -l | tr -d ' ')" 15
same "PERMANOVA residual df: plot 12, subplot 24" \
     "$(col "$P" residual_df distance bray_curtis term nutrient) $(col "$P" residual_df distance bray_curtis term exclosure)" "12 24"
same "plot-stratum units are 18 plot means" "$(col "$P" n_units distance bray_curtis term nutrient)" 18
within "nutrient effect found among plots (Bray-Curtis p)" "$(col "$P" p distance bray_curtis term nutrient)" 0 0.05
within "exclosure effect found within plots (Bray-Curtis p)" "$(col "$P" p distance bray_curtis term exclosure)" 0 0.05
within "depth effect found within subplots (Bray-Curtis p)" "$(col "$P" p distance bray_curtis term depth_layer)" 0 0.05
within "no round effect (Bray-Curtis p)" "$(col "$P" p distance bray_curtis term round)" 0.05 1
same "interaction term kept under its design name" "$(col "$P" stratum distance jaccard term nutrient:exclosure)" subplot
C="$A/da_consensus.tsv"
same "DA consensus: full exclosure -> Genus04-06" "$(genera "$C" exclosurefull_exclosure)" "Genus04 Genus05 Genus06"
same "DA consensus: 7 cm -> Genus07-09" "$(genera "$C" depth_layer7cm)" "Genus07 Genus08 Genus09"
same "DA consensus: no taxon outside the true effects" \
     "$(col "$C" feature | sed 's/.*g__//' | grep -cv '^Genus0[1-9]$' || true)" 0
same "coefficient meaning spells out the simple effect" \
     "$(col "$A/da_coefficients.tsv" meaning coefficient nutrientenriched)" \
     "nutrient: enriched vs unenriched_control, at exclosure = open_control"
same "ALDEx2 kept as supporting evidence" "$(col "$A/da_all_methods.tsv" role method ALDEx2 | sort -u)" supporting

# depth rules agree with pipeline/choose_depth.py on the same samples
python3 - "$WORK/data/feature-table.tsv" "$A/samples_removed.tsv" > "$WORK/depths.tsv" <<'EOF'
import csv, sys
rows = list(csv.reader(open(sys.argv[1]), delimiter='\t'))
removed = {r[0] for r in list(csv.reader(open(sys.argv[2]), delimiter='\t'))[1:]}
names = rows[1][1:]
tot = [sum(int(r[i + 1]) for r in rows[2:]) for i in range(len(names))]
print('sample\treads')
for n, t in zip(names, tot):
    if n not in removed:
        print(f'{n}\t{t}')
EOF
py="$(python3 "$REPO/pipeline/choose_depth.py" --depths "$WORK/depths.tsv" --metadata "$WORK/data/metadata.tsv" \
      --groups "round depth_layer" --out "$WORK/dc.tsv")"
same "R and choose_depth.py pick the same main and sensitivity depths" \
     "$(col "$A/depth_choice.tsv" value setting main_depth) $(col "$A/depth_choice.tsv" value setting sensitivity_depth)" "$py"
same "sensitivity re-run written" \
     "$([[ -s "$A/alpha_lmm_sensitivity.tsv" && -s "$A/permanova_hierarchical_sensitivity.tsv" ]] && echo yes)" yes

# ---- B. Manuscript-1 design: 0 cm only -> one sample per subplot ----
design "$WORK/ms1.tsv" ms1_photosynthetic_0cm.tsv
Rscript "$REPO/R/design_analysis.R" "$WORK/ms1.tsv" "$WORK/ms1" > "$WORK/ms1.log" 2>&1 \
    || { echo "FAIL  design_analysis.R (ms1) failed: $(tail -5 "$WORK/ms1.log")"; exit 1; }
B="$WORK/ms1"
same "0 cm filter keeps 54 samples" "$(grep -c "^GA.\{3\}0S" <(col "$B/alpha_diversity.tsv" sample))" 54
same "one sample per subplot: exclosure tested on samples within plots" \
     "$(col "$B/design_strata.tsv" tested_on term exclosure)" "samples (one per subplot)"
same "subplot random intercept dropped (confounded with residual)" \
     "$(col "$B/design_used.tsv" value setting random_effects)" "(1 | plot)"
within "alpha LMM: exclosure df at 0 cm" "$(col "$B/alpha_lmm.tsv" DenDF metric Shannon term exclosure)" 22 26
same "PERMANOVA residual df at 0 cm: plot 12, subplot 24" \
     "$(col "$B/permanova_hierarchical.tsv" residual_df distance bray_curtis term nutrient) $(col "$B/permanova_hierarchical.tsv" residual_df distance bray_curtis term exclosure)" "12 24"

# ---- C. A design file whose expected strata do not match the metadata stops ----
design "$WORK/bad.tsv" ms1_photosynthetic_0cm.tsv -e 's/^expect_strata\t.*/expect_strata\tnutrient=subplot/'
rc=0
Rscript "$REPO/R/design_analysis.R" "$WORK/bad.tsv" "$WORK/bad" > "$WORK/bad.log" 2>&1 || rc=$?
same "expect_strata mismatch stops with a message" \
     "$([[ $rc -ne 0 ]] && echo stopped) $(grep -c "varies at the plot level, but the design file expects subplot" "$WORK/bad.log")" \
     "stopped 1"

exit "$status"
