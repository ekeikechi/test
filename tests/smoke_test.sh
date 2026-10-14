#!/usr/bin/env bash
# Dry-run every pipeline branch against stub qiime/biom commands (no QIIME 2 needed).
# Checks control flow, file chaining between stages, checkpoints, the preset
# mechanism and the per-stage settings report. It does NOT test any results.
#
# Usage: bash tests/smoke_test.sh WORKDIR     (WORKDIR must not exist yet)
# To test an older bash (e.g. 4.3 on an HPC), put it first on PATH as `bash`.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/smoke_test.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
mkdir -p "$WORK"
export PATH="${REPO}/tests/stubs:${PATH}"
echo "bash under test: $(bash --version | head -n1)"
status=0

# ---- 1. Preset block: each "shipped" value must equal the value assigned above it ----
sed 's/^PRESET=""/PRESET="18S_PR2"/' "${REPO}/pipeline/config.sh" > "$WORK/config_check.sh"
mismatches="$(bash -c '
    set -u
    preset() { [[ "${!1-}" == "$2" ]] || echo "  preset $1: shipped value \"$2\" but config.sh has \"${!1-}\""; }
    override() { :; }
    source "$1"' _ "$WORK/config_check.sh")"
if [[ -z "$mismatches" ]]; then echo "PASS  preset block matches the shipped config values"
else echo "FAIL  preset block out of sync with config.sh:"; echo "$mismatches"; status=1; fi

# ---- 2. Branches: name|config edits. VAR=value replaces that line in config.sh;
#      an item starting with "+" is appended as a line (the "Your overrides" section). ----
BRANCHES=(
  'illumina16S_2026.7'
  'illumina16S_2026.4|Q2_EPOCH="2026.4"'
  'deblur18S_2026.1|DENOISER="deblur"|MARKER="18S"|Q2_EPOCH="2026.1"'
  'pacbioITS|PLATFORM="pacbio"|MARKER="ITS"|IMPORT_MODE="manifest"|RUN_KMER_DIVERSITY="yes"'
  'illumina_single|LAYOUT="single"|IMPORT_MODE="manifest"|PRIMER_R=""'
  'no_primers|PRIMER_F=""|PRIMER_R=""'
  'pr2_18S_preset|PRESET="18S_PR2"|CLASSIFIER=""'
  'pr2_user_override|PRESET="18S_PR2"|CLASSIFIER=""|MIN_SAMPLES=3|TAXA_COLLAPSE_LEVEL=7'
  'pr2_override_section|PRESET="18S_PR2"|CLASSIFIER=""|+override MIN_SAMPLES "2"'
  'pr2_auto_trunc|PRESET="18S_PR2"|CLASSIFIER=""|DADA2_TRUNC_F=""|DADA2_TRUNC_R=""'
  'auto_trunc_lengthen|DADA2_TRUNC_F="auto"|DADA2_TRUNC_R="auto"|DADA2_TRUNC_SUM_TARGET=520|DADA2_MAX_EE_R=3'
  'depth_auto|MIN_SAMPLE_DEPTH=1000|SAMPLING_DEPTH="auto"|DEPTH_SENSITIVITY_GROUPS="grp"|DEPTH_SENSITIVITY_MIN_KEEP=0.5'
)
# quality profile the stub exports for 02_trim/reads_trimmed_summary.qzv (tests/make_seven_number.py)
profile_for() { case "$1" in auto_trunc_lengthen) echo early_drop ;; *) echo 2x300 ;; esac; }
STAGES=(01_import_qc 02_trim_primers 03_denoise 04_filter_features 05a_build_pr2_classifier
        05_taxonomy 06_filter_taxa_samples 07_phylogeny 08_diversity 09_composition_da
        10_export_provenance)

# expect NAME DESCRIPTION COMMAND...: an assertion for a branch
expect() {
    local name="$1" what="$2"; shift 2
    if "$@" >/dev/null 2>&1; then echo "      ok   $what"; else echo "      FAIL $what"; status=1; fi
}
# shellcheck disable=SC2329  # called through expect "$@"
setting_is() {  # setting_is STAGE_DIR VAR VALUE SOURCE_PREFIX
    awk -F'\t' -v v="$2" -v val="$3" -v s="$4" \
        '$1 == v && $2 == val && index($3, s) == 1 { found = 1 } END { exit !found }' \
        "$B/results/$1/settings_used.tsv"
}

for spec in "${BRANCHES[@]}"; do
    IFS='|' read -r -a parts <<< "$spec"
    name="${parts[0]}"; B="${WORK}/${name}"
    mkdir -p "$B/raw_data"
    cp -r "${REPO}/pipeline" "$B/"
    touch "$B/manifest.tsv" "$B/raw_data/a_S1_L001_R1_001.fastq.gz"
    # metadata with a group column, except one branch with IDs only (as at the start of a project)
    # (samples a-d match the stub biom's per-sample counts: 1200, 2500, 4000, 900)
    if [[ "$name" == pr2_auto_trunc ]]; then printf 'sample-id\na\nb\nc\nd\n' > "$B/metadata.tsv"
    else printf 'sample-id\tgrp\na\tx\nb\tx\nc\ty\nd\ty\n' > "$B/metadata.tsv"; fi
    for f in clf ref sepp; do echo stub > "$B/$f.qza"; done
    cfg="$B/pipeline/config.sh"
    # answer every checkpoint, then apply the branch's own settings
    for kv in "CLASSIFIER=\"$B/clf.qza\"" "SEPP_REF=\"$B/sepp.qza\"" "DEBLUR_REF_SEQS=\"$B/ref.qza\"" \
              'MIN_SAMPLE_DEPTH=1' 'SAMPLING_DEPTH=1000' 'GROUP_COL="grp"' 'ADONIS_FORMULA="grp"' \
              'DADA2_TRUNC_F=150' 'DADA2_TRUNC_R=150' 'DEBLUR_TRIM_LENGTH=250' "${parts[@]:1}"; do
        if [[ "$kv" == +* ]]; then echo "${kv#+}" >> "$cfg"
        else sed -i "s|^${kv%%=*}=.*|${kv}|" "$cfg"; fi
    done
    export SMOKE_LOG="$B/calls.log"; : > "$SMOKE_LOG"
    STUB_EPOCH="$(sed -n 's/^Q2_EPOCH="\(.*\)".*/\1/p' "$cfg")"; export STUB_EPOCH
    STUB_QUALITY_PROFILE="$(profile_for "$name")"; export STUB_QUALITY_PROFILE
    result="OK"
    for s in "${STAGES[@]}"; do
        if ! (cd "$B" && bash "pipeline/${s}.sh") > "$B/${s}.log" 2>&1; then
            result="FAILED at ${s} (see $B/${s}.log)"; status=1; break
        fi
    done
    printf '%-6s%-20s (%s qiime calls)\n' "$([[ $result == OK ]] && echo PASS || echo FAIL)" \
        "$name" "$(grep -c '^qiime' "$SMOKE_LOG")"
    [[ "$result" == OK ]] || { echo "      $result"; continue; }

    expect "$name" "every stage printed and saved its settings" \
        test "$(find "$B/results" -name settings_used.tsv | wc -l)" -eq 11
    case "$name" in
        illumina16S_2026.7)
            expect "$name" "collapse at level 6 (genus)" grep -q -- '--p-level 6 --o-collapsed-table .*table_da_genus.qza' "$SMOKE_LOG"
            expect "$name" "05a skipped (BUILD_PR2_CLASSIFIER=no)" bash -c "! grep -q 'rescript' '$SMOKE_LOG'"
            expect "$name" "settings handed to R" grep -qx $'collapse_level\t6' "$B/results/10_export/analysis_settings.tsv"
            expect "$name" "no PR2 version without a PR2 classifier" grep -qx $'pr2_version\t' "$B/results/10_export/analysis_settings.tsv"
            expect "$name" "metadata with columns: summarize and grouped curves use it" \
                bash -c "grep -q 'feature-table summarize .*--m-metadata-file' '$SMOKE_LOG' && grep -q 'alpha-rarefaction .*--m-metadata-file' '$SMOKE_LOG'"
            expect "$name" "cutadapt --o-stats, tabulated and exported (2026.7)" \
                bash -c "grep -q -- '--o-stats .*02_trim/cutadapt_stats.qza' '$SMOKE_LOG' && grep -q -- 'tools export --input-path .*cutadapt_stats.qza' '$SMOKE_LOG'"
            expect "$name" "da-barplot for ANCOM-BC2 (2026.7)" grep -q 'composition da-barplot' "$SMOKE_LOG" ;;
        illumina16S_2026.4)
            expect "$name" "no --o-stats before 2026.7" bash -c "! grep -q -- 'cutadapt_stats' '$SMOKE_LOG'"
            expect "$name" "ancombc2-visualizer on 2026.4" bash -c "grep -q 'ancombc2-visualizer' '$SMOKE_LOG' && ! grep -q 'da-barplot' '$SMOKE_LOG'" ;;
        pr2_18S_preset)
            expect "$name" "PR2 classifier built (get-pr2-data, dereplicate --p-rank-handles disable, fit)" \
                bash -c "grep -q 'rescript get-pr2-data --p-version 5.1.1' '$SMOKE_LOG' && grep -q -- '--p-rank-handles disable' '$SMOKE_LOG' && grep -q 'fit-classifier-naive-bayes' '$SMOKE_LOG'"
            expect "$name" "stage 05 used the built PR2 classifier" grep -q -- '--i-classifier .*05a_pr2_classifier/pr2_5.1.1_classifier.qza' "$SMOKE_LOG"
            expect "$name" "18S V4 primers" grep -q -- '--p-front-f CYGCGGTAATTCCAGCTC --p-front-r AYGGTATCTRATCRTCTTYG' "$SMOKE_LOG"
            expect "$name" "DADA2 pseudo-pooling" grep -q -- '--p-pooling-method pseudo' "$SMOKE_LOG"
            expect "$name" "MIN_SAMPLES=1 from preset" setting_is 04_filter_features MIN_SAMPLES 1 "preset 18S_PR2"
            expect "$name" "MIN_FREQ_FRACTION=0.001 (keep 0.1% rule) from preset" setting_is 04_filter_features MIN_FREQ_FRACTION 0.001 "preset 18S_PR2"
            expect "$name" "PR2 exclusions on domain + organelle tags" grep -q -- '--p-exclude d__Bacteria,d__Archaea,:plas,:mito,:apic,:chro' "$SMOKE_LOG"
            expect "$name" "no tree; non-phylogenetic diversity" bash -c "! grep -q -- '--i-phylogeny\|fragment-insertion\|align-to-tree' '$SMOKE_LOG'"
            expect "$name" "kmer diversity on" grep -q 'boots kmer-diversity' "$SMOKE_LOG"
            expect "$name" "collapse at level 8 (genus of 9 PR2 ranks)" grep -q -- '--p-level 8 --o-collapsed-table .*table_da_genus.qza' "$SMOKE_LOG"
            expect "$name" "R/pr2_traits.R gets the classifier's PR2 version" \
                grep -qx $'pr2_version\t5.1.1' "$B/results/10_export/analysis_settings.tsv"
            expect "$name" "R gets 9 ranks + level 8" bash -c "grep -qx \$'taxonomy_ranks\tDomain Supergroup Division Subdivision Class Order Family Genus Species' '$B/results/10_export/analysis_settings.tsv' && grep -qx \$'collapse_level\t8' '$B/results/10_export/analysis_settings.tsv'" ;;
        pr2_user_override)
            expect "$name" "your MIN_SAMPLES=3 kept and labelled" setting_is 04_filter_features MIN_SAMPLES 3 "yours; preset 18S_PR2 would set '1'"
            expect "$name" "your TAXA_COLLAPSE_LEVEL=7 kept (family)" grep -q -- '--p-level 7 --o-collapsed-table .*table_da_family.qza' "$SMOKE_LOG"
            expect "$name" "other preset values still applied" setting_is 06_filter_taxa_samples TAXA_EXCLUDE "d__Bacteria,d__Archaea,:plas,:mito,:apic,:chro" "preset 18S_PR2" ;;
        pr2_override_section)
            expect "$name" "override keeps MIN_SAMPLES=2 although it equals the 16S default" \
                setting_is 04_filter_features MIN_SAMPLES 2 "yours: override; preset 18S_PR2 would set '1'"
            expect "$name" "stage 04 filtered with --p-min-samples 2" grep -q -- '--p-min-samples 2 ' "$SMOKE_LOG"
            expect "$name" "rest of the preset still applied" setting_is 04_filter_features MIN_FREQ_FRACTION 0.001 "preset 18S_PR2" ;;
        pr2_auto_trunc)
            expect "$name" "preset sets auto; quality rule gives 282/250 (sum 532 >= 520)" \
                grep -q -- '--p-trunc-len-f 282 --p-trunc-len-r 250 ' "$SMOKE_LOG"
            expect "$name" "settings report shows the chosen lengths, marked auto" \
                setting_is 03_denoise DADA2_TRUNC_R 250 "auto, chosen in truncation_auto.tsv; preset 18S_PR2"
            expect "$name" "decision table written" grep -qx $'status\tok' "$B/results/03_denoise/truncation_auto.tsv"
            expect "$name" "ID-only metadata: summarize runs without it" \
                bash -c "grep 'feature-table summarize' '$SMOKE_LOG' | grep -vq -- '--m-metadata-file'"
            expect "$name" "ID-only metadata: grouped rarefaction curves skipped" \
                bash -c "[[ \$(grep -c 'diversity alpha-rarefaction' '$SMOKE_LOG') == 1 ]] && ! grep -q 'alpha-rarefaction .*--m-metadata-file' '$SMOKE_LOG'"
            expect "$name" "18S_PR2 preset: expected-error filter 2/3" grep -q -- '--p-max-ee-f 2 --p-max-ee-r 3 ' "$SMOKE_LOG"
            expect "$name" "post-DADA2 check flags the bad sample" \
                bash -c "grep -q \$'^S3\t.*low_filter_pass' '$B/results/03_denoise/denoising_flags.tsv' && ! grep -q \$'^S1\t.*low' '$B/results/03_denoise/denoising_flags.tsv'" ;;
        auto_trunc_lengthen)
            expect "$name" "reverse lengthened from 200 to 238 to reach F+R = 520" \
                grep -q -- '--p-trunc-len-f 282 --p-trunc-len-r 238 ' "$SMOKE_LOG"
            expect "$name" "DADA2_MAX_EE_R=3 passed to DADA2" grep -q -- '--p-max-ee-r 3 ' "$SMOKE_LOG"
            expect "$name" "lengthening recorded" grep -qx $'reverse_reason\tlengthened for overlap' "$B/results/03_denoise/truncation_auto.tsv" ;;
        depth_auto)
            expect "$name" "sample below MIN_SAMPLE_DEPTH listed as removed" grep -q $'^d\t900' "$B/results/06_filter_taxa_samples/samples_removed.tsv"
            expect "$name" "auto depth = smallest sample (900)" grep -q -- 'boots core-metrics .*--p-sampling-depth 900 ' "$SMOKE_LOG"
            expect "$name" "sensitivity run at 2500 (every group keeps >= 50%)" grep -q -- 'boots core-metrics .*--p-sampling-depth 2500 .*sensitivity/core' "$SMOKE_LOG"
            expect "$name" "per-group table lists what the sensitivity depth drops" grep -q $'^x\t2\t2\t1\ta' "$B/results/08_diversity/depth_by_group.tsv"
            expect "$name" "settings report shows the chosen depth" setting_is 08_diversity SAMPLING_DEPTH 900 "auto, chosen in depth_choice.tsv"
            expect "$name" "depth rules handed to R" grep -qx $'depth_sensitivity_groups\tgrp' "$B/results/10_export/analysis_settings.tsv" ;;
        illumina_single)
            expect "$name" "single-end gets --p-max-ee" grep -q -- 'denoise-single .*--p-max-ee 2 ' "$SMOKE_LOG" ;;
        deblur18S_2026.1)
            expect "$name" "2026.1 visualizer, no seed flag" bash -c "grep -q 'ancombc2-visualizer' '$SMOKE_LOG' && ! grep -q random-seed '$SMOKE_LOG'" ;;
        pacbioITS)
            expect "$name" "no tree for ITS; kmer diversity" bash -c "! grep -q -- '--i-phylogeny' '$SMOKE_LOG' && grep -q 'boots kmer-diversity' '$SMOKE_LOG'" ;;
    esac
done

# ---- 3. Configuration errors are caught before anything runs ----
B="${WORK}/bad_level"; mkdir -p "$B"; cp -r "${REPO}/pipeline" "$B/"
sed -i 's|^TAXA_COLLAPSE_LEVEL=.*|TAXA_COLLAPSE_LEVEL=8|' "$B/pipeline/config.sh"   # 8 > 7 ranks
if (cd "$B" && bash pipeline/01_import_qc.sh) > "$B/out.log" 2>&1; then
    echo "FAIL  collapse level beyond TAXONOMY_RANKS was not rejected"; status=1
elif grep -q 'TAXA_COLLAPSE_LEVEL=8 but TAXONOMY_RANKS lists only 7 ranks' "$B/out.log"; then
    echo "PASS  collapse level beyond TAXONOMY_RANKS rejected with a clear message"
else
    echo "FAIL  unexpected error: $(tail -2 "$B/out.log")"; status=1
fi

# ---- 4. Automatic truncation stops when the overlap target cannot be reached ----
B="${WORK}/auto_unreachable"; mkdir -p "$B/raw_data"; cp -r "${REPO}/pipeline" "$B/"
printf 'sample-id\na\n' > "$B/metadata.tsv"; touch "$B/raw_data/a_S1_L001_R1_001.fastq.gz"
sed -i 's|^PRESET=.*|PRESET="18S_PR2"|' "$B/pipeline/config.sh"
export SMOKE_LOG="$B/calls.log" STUB_EPOCH=2026.7 STUB_QUALITY_PROFILE=2x250; : > "$SMOKE_LOG"
if (cd "$B" && for s in 01_import_qc 02_trim_primers 03_denoise; do bash "pipeline/${s}.sh" || exit 1; done) > "$B/out.log" 2>&1; then
    echo "FAIL  2x250 reads with the 520 target were not stopped"; status=1
elif grep -q 'Overlap target not reachable: forward + reverse can reach at most 463 nt' "$B/out.log" \
     && ! grep -q 'dada2 denoise-paired' "$SMOKE_LOG"; then
    echo "PASS  unreachable overlap target stops stage 03 before DADA2, with a clear message"
else
    echo "FAIL  unexpected result: $(tail -3 "$B/out.log")"; status=1
fi

exit "$status"
