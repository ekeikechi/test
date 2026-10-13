# shellcheck shell=bash
# Shared helpers for the stage scripts. Sourced, never executed directly.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${PIPELINE_DIR}/config.sh"

log() { printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

# need VAR "message": stop with guidance when a decision variable is still empty.
need() {
    local name="$1" hint="$2"
    [[ -n "${!name:-}" ]] || die "${name} is not set in pipeline/config.sh. ${hint}"
}

# have FILE: true when a stage output already exists (lets stages be re-run safely).
have() { [[ -e "$1" ]]; }

# run CMD...: echo the command, then execute it.
run() { printf '+ %s\n' "$*" >&2; "$@"; }

# epoch_ge X: true when Q2_EPOCH >= X (both as YYYY.M).
epoch_ge() {
    local a="${Q2_EPOCH}" b="$1"
    [[ "$(printf '%s\n%s\n' "$a" "$b" | sort -V | head -n1)" == "$b" ]]
}

# depth_stat TABLE.qza EXPORT_DIR STAT: per-sample depth statistic (Min, Max, Median or
# Mean) from the biom CLI shipped with QIIME 2. The BIOM export is kept in EXPORT_DIR.
depth_stat() {
    local table="$1" dir="$2" stat="$3"
    have "$dir/feature-table.biom" \
        || qiime tools export --input-path "$table" --output-path "$dir" >/dev/null
    biom summarize-table -i "$dir/feature-table.biom" \
        | awk -F': ' -v s=" ${stat}" '$1 == s { gsub(",", "", $2); print $2; exit }'
}

# summarize TABLE.qza PREFIX: feature-table summary (2026.x pipeline interface).
summarize() {
    local table="$1" prefix="$2"
    have "${prefix}_summary.qzv" && return 0
    run qiime feature-table summarize \
        --i-table "$table" \
        --m-metadata-file "$METADATA" \
        --o-summary "${prefix}_summary.qzv" \
        --o-sample-frequencies "${prefix}_sample_frequencies.qza" \
        --o-feature-frequencies "${prefix}_feature_frequencies.qza"
}

command -v qiime >/dev/null 2>&1 || die "qiime not found: activate the QIIME 2 ${Q2_EPOCH} environment first."
mkdir -p "$OUT"
