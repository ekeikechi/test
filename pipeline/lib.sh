# shellcheck shell=bash
# Shared helpers for the stage scripts. Sourced, never executed directly.
#
# Portability: works on bash 3.2+. Arrays that may be empty are expanded as
# ${arr[@]+"${arr[@]}"}, because bash < 4.4 treats "${arr[@]}" of an empty array as
# an unbound variable under `set -u`.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# preset VAR SHIPPED VALUE: used by the preset block in config.sh.
# Sets VAR to VALUE only while VAR still holds the value shipped in config.sh, so a
# preset never overwrites a value you changed yourself. The outcome is recorded in
# PRESET_SRC_<VAR> and reported by show_settings.
preset() {
    local var="$1" shipped="$2" value="$3"
    if [[ "${!var-}" == "$shipped" ]]; then
        printf -v "$var" '%s' "$value"
        printf -v "PRESET_SRC_${var}" '%s' "preset ${PRESET}"
    else
        printf -v "PRESET_SRC_${var}" '%s' "yours; preset ${PRESET} would set '${value}'"
    fi
}

# override VAR VALUE: used in the "Your overrides" section at the end of config.sh.
# Applied after the preset, so it always wins; reported as (yours: override ...).
override() {
    local var="$1" value="$2" src_var="PRESET_SRC_$1" note="yours: override"
    case "${!src_var:-}" in
        "preset "*) note="yours: override; ${!src_var} would set '${!var-}'" ;;
    esac
    printf -v "$var" '%s' "$value"
    printf -v "PRESET_SRC_${var}" '%s' "$note"
}

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

# show_settings STAGE_DIR VAR...: print the settings this stage uses, each with its
# source (config / preset / yours), at the start of the log, and save them to
# STAGE_DIR/settings_used.tsv.
show_settings() {
    local dir="$1" var src_var
    shift
    mkdir -p "$dir"
    {
        printf 'setting\tvalue\tsource\n'
        for var in "$@"; do
            src_var="PRESET_SRC_${var}"
            printf '%s\t%s\t%s\n' "$var" "${!var-}" "${!src_var:-config}"
        done
    } > "$dir/settings_used.tsv"
    printf '\n[%s] Settings used by %s (preset: %s)\n' \
        "$(date '+%H:%M:%S')" "$(basename "$0")" "${PRESET:-none}" >&2
    awk -F'\t' 'NR > 1 { printf "  %-30s = %-40s (%s)\n", $1, ($2 == "" ? "<empty>" : $2), $3 }' \
        "$dir/settings_used.tsv" >&2
}

# rank_count / rank_name N: the configured taxonomy ranks (TAXONOMY_RANKS, 1-based).
rank_count() {
    # shellcheck disable=SC2086  # word splitting of the rank list is intended
    set -- $TAXONOMY_RANKS
    echo "$#"
}
rank_name() {
    local i="$1"
    # shellcheck disable=SC2086
    set -- $TAXONOMY_RANKS
    echo "${!i}"
}

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

# metadata_columns: number of data columns (besides the sample ID) in $METADATA.
# Header = first non-blank line that is not a comment ("#SampleID"-style headers count).
metadata_columns() {
    awk -F'\t' '/^[[:space:]]*$/ { next }
        /^#/ && $1 !~ /^#(SampleID|Sample ID|OTUID|OTU ID|sample_name)$/ { next }
        { n = NF - 1; exit } END { print n + 0 }' "$METADATA"
}

# sample_depths TABLE.qza EXPORT_DIR OUT.tsv: reads per sample ("sample<TAB>reads"), from the
# "Counts/sample detail" list of biom summarize-table. The BIOM export is kept in EXPORT_DIR.
sample_depths() {
    local table="$1" dir="$2" out="$3"
    have "$dir/feature-table.biom" \
        || qiime tools export --input-path "$table" --output-path "$dir" >/dev/null
    biom summarize-table -i "$dir/feature-table.biom" | awk '
        BEGIN { print "sample\treads" }
        /^Counts\/sample detail:/ { on = 1; next }
        on && /:/ { i = index($0, ":"); v = substr($0, i + 1); gsub(/[ ,]/, "", v)
                    print substr($0, 1, i - 1) "\t" v }' > "$out"
}

# summarize TABLE.qza PREFIX: feature-table summary (2026.x pipeline interface).
# Metadata is passed only when it has a data column: summarize (2026.7) fails with
# "list index out of range" on ID-only metadata.
summarize() {
    local table="$1" prefix="$2" meta_args=()
    have "${prefix}_summary.qzv" && return 0
    (( $(metadata_columns) > 0 )) && meta_args=(--m-metadata-file "$METADATA")
    run qiime feature-table summarize \
        --i-table "$table" \
        ${meta_args[@]+"${meta_args[@]}"} \
        --o-summary "${prefix}_summary.qzv" \
        --o-sample-frequencies "${prefix}_sample_frequencies.qza" \
        --o-feature-frequencies "${prefix}_feature_frequencies.qza"
}

# ---- configuration checks ----------------------------------------------------
case "${PRESET}" in ""|18S_PR2) ;; *) die "Unknown PRESET '${PRESET}' (use \"\" or \"18S_PR2\")." ;; esac
for level_var in TAXA_COLLAPSE_LEVEL OVERVIEW_RANK_LEVEL; do
    level="${!level_var}"
    [[ "$level" =~ ^[1-9][0-9]*$ ]] || die "${level_var} must be a positive integer (got '${level}')."
    (( level <= $(rank_count) )) \
        || die "${level_var}=${level} but TAXONOMY_RANKS lists only $(rank_count) ranks: ${TAXONOMY_RANKS}"
done
unset level_var level

command -v qiime >/dev/null 2>&1 || die "qiime not found: activate the QIIME 2 ${Q2_EPOCH} environment first."
mkdir -p "$OUT"
