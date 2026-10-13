#!/usr/bin/env bash
# Dry-run every pipeline branch against stub qiime/biom commands (no QIIME 2 needed).
# Usage: bash tests/smoke_test.sh WORKDIR     (WORKDIR must not exist yet)
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:?usage: bash tests/smoke_test.sh WORKDIR}"
[[ -e "$WORK" ]] && { echo "WORKDIR already exists: $WORK" >&2; exit 1; }
export PATH="${REPO}/tests/stubs:${PATH}"

# name followed by config overrides (VAR=value, written verbatim into config.sh)
BRANCHES=(
  'illumina16S_2026.7'
  'deblur18S_2026.1|DENOISER="deblur"|MARKER="18S"|Q2_EPOCH="2026.1"'
  'pacbioITS|PLATFORM="pacbio"|MARKER="ITS"|IMPORT_MODE="manifest"|RUN_KMER_DIVERSITY="yes"'
  'illumina_single|LAYOUT="single"|IMPORT_MODE="manifest"|PRIMER_R=""'
  'no_primers|PRIMER_F=""|PRIMER_R=""'
)
STAGES=(01_import_qc 02_trim_primers 03_denoise 04_filter_features 05_taxonomy
        06_filter_taxa_samples 07_phylogeny 08_diversity 09_composition_da 10_export_provenance)

status=0
for spec in "${BRANCHES[@]}"; do
    IFS='|' read -r -a parts <<< "$spec"
    name="${parts[0]}"; B="${WORK}/${name}"
    mkdir -p "$B/raw_data"
    cp -r "${REPO}/pipeline" "$B/"
    touch "$B/metadata.tsv" "$B/manifest.tsv" "$B/raw_data/a_S1_L001_R1_001.fastq.gz"
    for f in clf ref sepp; do echo stub > "$B/$f.qza"; done
    cfg="$B/pipeline/config.sh"
    # answer every checkpoint
    for kv in "CLASSIFIER=\"$B/clf.qza\"" "SEPP_REF=\"$B/sepp.qza\"" "DEBLUR_REF_SEQS=\"$B/ref.qza\"" \
              'MIN_SAMPLE_DEPTH=1' 'SAMPLING_DEPTH=1000' 'GROUP_COL="grp"' 'ADONIS_FORMULA="grp"' \
              'DADA2_TRUNC_F=150' 'DADA2_TRUNC_R=150' 'DEBLUR_TRIM_LENGTH=250' "${parts[@]:1}"; do
        sed -i "s|^${kv%%=*}=.*|${kv}|" "$cfg"
    done
    export SMOKE_LOG="$B/calls.log"; : > "$SMOKE_LOG"
    result="OK"
    for s in "${STAGES[@]}"; do
        if ! (cd "$B" && bash "pipeline/${s}.sh") > "$B/${s}.log" 2>&1; then
            result="FAILED at ${s} (see $B/${s}.log)"; status=1; break
        fi
    done
    printf '%-20s %s  (%s qiime calls)\n' "$name" "$result" "$(grep -c '^qiime' "$SMOKE_LOG")"
done
exit "$status"
