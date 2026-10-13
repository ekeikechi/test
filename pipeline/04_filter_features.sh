#!/usr/bin/env bash
# Stage 04: remove rare ASVs before taxonomy (faster, less memory).
# Rules combined from:
#   prevalence  >= MIN_SAMPLES samples           [gut_to_soil_tutorial: filter_features min_samples=2]
#   frequency   >= MIN_FREQ_FRACTION x mean depth [MH_qiime2-amplicon-2024.5 §4.1] [MH_2026.3-workflow §4.1]
#   filter before classifying                     [gut_to_soil_tutorial] [conceptual_overview]
source "$(dirname "$0")/lib.sh"
D="${OUT}/04_filter_features"; mkdir -p "$D"
show_settings "$D" \
    MIN_SAMPLES MIN_FREQ_FRACTION
IN_TABLE="${OUT}/03_denoise/table.qza"
IN_SEQS="${OUT}/03_denoise/rep_seqs.qza"

MEAN_DEPTH="$(depth_stat "$IN_TABLE" "$D/denoised_table_biom" Mean)"
# round(fraction x mean depth), done in awk to avoid a bc dependency.
MIN_FREQ="$(awk -v m="$MEAN_DEPTH" -v f="$MIN_FREQ_FRACTION" 'BEGIN { printf "%d", m * f + 0.5 }')"
log "Mean sample depth = ${MEAN_DEPTH}; min ASV frequency = ${MIN_FREQ}; min samples = ${MIN_SAMPLES}"
printf 'mean_depth\t%s\nmin_frequency\t%s\nmin_samples\t%s\n' \
    "$MEAN_DEPTH" "$MIN_FREQ" "$MIN_SAMPLES" > "$D/filter_thresholds.tsv"

have "$D/table.qza" || run qiime feature-table filter-features \
    --i-table "$IN_TABLE" \
    --p-min-samples "$MIN_SAMPLES" \
    --p-min-frequency "$MIN_FREQ" \
    --o-filtered-table "$D/table.qza"

have "$D/rep_seqs.qza" || run qiime feature-table filter-seqs \
    --i-data "$IN_SEQS" \
    --i-table "$D/table.qza" \
    --o-filtered-data "$D/rep_seqs.qza"

summarize "$D/table.qza" "$D/table"

log "CHECKPOINT: compare $D/table_summary.qzv with 03_denoise/table_summary.qzv."
log "Expect many ASVs but few reads to be removed. If a large share of reads is lost,"
log "relax MIN_SAMPLES (e.g. studies with unreplicated, unique samples) and re-run."
