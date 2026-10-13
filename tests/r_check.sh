#!/usr/bin/env bash
# Run R/downstream_analysis.R on synthetic data with 6 known enriched genera and check
# the results. Scenarios:
#   A  16S 7-rank, no settings file (script defaults)           -> consensus = Genus01-06
#   B  16S, DA_REFERENCE "treatment::treated" (non-alphabetical) -> consensus = Genus01-06
#      (regression test: ALDEx2's sign must be aligned with the reference)
#   C  PR2 9-rank, settings from the 18S_PR2 preset (level 8)    -> consensus = Genus01-06,
#      DA run at "Genus", empty multi-letter labels (dvs__) read as unassigned
#   D  PR2 taxonomy but 7 configured ranks                       -> R must stop, not truncate
#   E  read_taxonomy unit check: "sgr__", "dvs__", "g__" all become NA
# Usage: bash tests/r_check.sh WORKDIR     (WORKDIR must not exist yet; needs Rscript + packages)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/r_check.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
EXPECTED="Genus01 Genus02 Genus03 Genus04 Genus05 Genus06 "
status=0

settings() {  # settings FILE RANKS COLLAPSE OVERVIEW REFERENCE
    printf 'setting\tvalue\npreset\t\ntaxonomy_ranks\t%s\ncollapse_level\t%s\noverview_level\t%s\nsampling_depth\t\ngroup_col\ttreatment\nda_reference\t%s\nrandom_seed\t711\n' \
        "$2" "$3" "$4" "$5" > "$1"
}
consensus() {
    awk -F, 'NR > 1 { gsub(/"/, "", $1); sub(/.*g__/, "", $1); print $1 }' "$1/da_consensus.csv" \
        | sort | tr '\n' ' '
}
report() {  # report NAME CONDITION_OK DETAIL
    if [[ "$2" == yes ]]; then echo "PASS  $1  $3"; else echo "FAIL  $1  $3"; status=1; fi
}
run_r() {  # run_r EXPORT OUT -> exit status of the R script
    Rscript "${REPO}/R/downstream_analysis.R" "$1" "$2" > "$2.log" 2>&1
}

R7="Domain Phylum Class Order Family Genus Species"
R9="Domain Supergroup Division Subdivision Class Order Family Genus Species"

# A
mkdir -p "$WORK/A"; Rscript "${REPO}/tests/make_synthetic.R" "$WORK/A" 16S >/dev/null
if run_r "$WORK/A" "$WORK/A_out"; then
    got="$(consensus "$WORK/A_out")"; [[ "$got" == "$EXPECTED" ]] && ok=yes || ok=no
    report "A 16S defaults" "$ok" "consensus: $got"
else report "A 16S defaults" no "R failed: $(tail -3 "$WORK/A_out.log")"; fi

# B
mkdir -p "$WORK/B"; Rscript "${REPO}/tests/make_synthetic.R" "$WORK/B" 16S >/dev/null
settings "$WORK/B/analysis_settings.tsv" "$R7" 6 2 "treatment::treated"
if run_r "$WORK/B" "$WORK/B_out"; then
    got="$(consensus "$WORK/B_out")"
    ref="$(awk -F, '$1 ~ /group_ref/ { gsub(/"/, "", $2); print $2 }' "$WORK/B_out/settings_used.csv")"
    [[ "$got" == "$EXPECTED" && "$ref" == "treated" ]] && ok=yes || ok=no
    report "B reference=treated" "$ok" "group_ref=$ref; consensus: $got"
else report "B reference=treated" no "R failed: $(tail -3 "$WORK/B_out.log")"; fi

# C
mkdir -p "$WORK/C"; Rscript "${REPO}/tests/make_synthetic.R" "$WORK/C" PR2 >/dev/null
settings "$WORK/C/analysis_settings.tsv" "$R9" 8 3 ""
if run_r "$WORK/C" "$WORK/C_out"; then
    got="$(consensus "$WORK/C_out")"
    rank="$(grep -o 'DA at [A-Za-z]* level' "$WORK/C_out.log" || true)"
    plots="$(cd "$WORK/C_out" && printf '%s ' barplot_*.pdf)"
    bogus="$(grep -c '"dvs__"' "$WORK/C_out/da_all_methods.csv" || true)"
    [[ "$got" == "$EXPECTED" && "$rank" == "DA at Genus level" && "$plots" == "barplot_division.pdf barplot_genus.pdf " && "$bogus" == 0 ]] && ok=yes || ok=no
    report "C PR2 9 ranks" "$ok" "${rank}; plots: ${plots}; consensus: $got"
else report "C PR2 9 ranks" no "R failed: $(tail -3 "$WORK/C_out.log")"; fi

# D
mkdir -p "$WORK/D"; Rscript "${REPO}/tests/make_synthetic.R" "$WORK/D" PR2 >/dev/null
settings "$WORK/D/analysis_settings.tsv" "$R7" 6 2 ""
if run_r "$WORK/D" "$WORK/D_out"; then
    report "D PR2 vs 7 ranks" no "R ran although the taxonomy is deeper than the ranks"
else
    msg="$(grep -o 'taxonomy.tsv has up to [0-9]* ranks but only [0-9]* rank names are configured' "$WORK/D_out.log" || true)"
    [[ -n "$msg" ]] && ok=yes || ok=no
    report "D PR2 vs 7 ranks" "$ok" "stopped: ${msg:-unexpected error}"
fi

# E: evaluate only the read_taxonomy() definition from the script
printf 'Feature ID\tTaxon\tConfidence\nf1\td__Eukaryota; sgr__; dv__X; dvs__; c__C1; o__; f__; g__; s__;\t1\n' > "$WORK/E_taxonomy.tsv"
e_out="$(Rscript -e '
  ex <- parse(commandArgs(TRUE)[1])
  for (e in ex) if (is.call(e) && length(e) == 3 && identical(e[[2]], as.name("read_taxonomy"))) eval(e)
  m <- read_taxonomy(commandArgs(TRUE)[2], c("Domain","Supergroup","Division","Subdivision",
                                            "Class","Order","Family","Genus","Species"))
  cat(paste(ifelse(is.na(m[1, ]), "NA", m[1, ]), collapse = " "))' \
  "${REPO}/R/downstream_analysis.R" "$WORK/E_taxonomy.tsv")"
[[ "$e_out" == "d__Eukaryota NA dv__X NA c__C1 NA NA NA NA" ]] && ok=yes || ok=no
report "E empty-label regex" "$ok" "parsed: $e_out"

exit "$status"
