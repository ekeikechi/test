#!/usr/bin/env bash
# Stage 06: remove non-target taxa, review depth, remove failed samples.
# Sources: [MH_qiime2-amplicon-2024.5 §4.2-4.4] [MH_2026.3-workflow §4.2-4.4]
#          [gut_to_soil_tutorial: remove failed samples; check per-group retention]
#          [MH_qiime2-amplicon-2024.5 §6: per-sample and per-group rarefaction curves]
source "$(dirname "$0")/lib.sh"
D="${OUT}/06_filter_taxa_samples"; mkdir -p "$D"
show_settings "$D" \
    MARKER TAXA_EXCLUDE TAXA_INCLUDE MIN_SAMPLE_DEPTH GROUP_COL
TABLE="${OUT}/04_filter_features/table.qza"
SEQS="${OUT}/04_filter_features/rep_seqs.qza"
TAX="${OUT}/05_taxonomy/taxonomy.qza"

# 6a. Taxonomy-based filter. Matching is case-insensitive substring ("contains").
# Unclassified ASVs are kept by default: many of them flag upstream problems or novelty.
if ! have "$D/table_taxfilt.qza"; then
    if [[ -n "$TAXA_EXCLUDE" || -n "$TAXA_INCLUDE" ]]; then
        tax_args=(--i-table "$TABLE" --i-taxonomy "$TAX" --o-filtered-table "$D/table_taxfilt.qza")
        [[ -n "$TAXA_EXCLUDE" ]] && tax_args+=(--p-exclude "$TAXA_EXCLUDE")
        [[ -n "$TAXA_INCLUDE" ]] && tax_args+=(--p-include "$TAXA_INCLUDE")
        run qiime taxa filter-table "${tax_args[@]}"
    else
        cp "$TABLE" "$D/table_taxfilt.qza"
    fi
fi
summarize "$D/table_taxfilt.qza" "$D/table_taxfilt"

# 6b. Rarefaction curves to pick a minimum depth (and later the even-sampling depth).
# Max depth = median sample depth works well [moving_pictures_tutorial] [gut_to_soil_tutorial];
# curves are drawn twice: grouped by metadata, and per sample (no metadata) [MH 2024.5 §6].
if ! have "$D/rarefaction_per_sample.qzv"; then
    MEDIAN="$(depth_stat "$D/table_taxfilt.qza" "$D/table_taxfilt_biom" Median)"
    MEDIAN="${MEDIAN%.*}"
    log "Rarefaction curves up to the median depth (${MEDIAN})"
    # Grouped curves need a categorical column; alpha-rarefaction fails on metadata
    # without one ("All metadata filtered after dropping columns ...").
    if [[ -n "$GROUP_COL" ]] && (( $(metadata_columns) > 0 )); then
        run qiime diversity alpha-rarefaction \
            --i-table "$D/table_taxfilt.qza" \
            --p-max-depth "$MEDIAN" --p-steps 20 \
            --p-metrics observed_features shannon \
            --m-metadata-file "$METADATA" \
            --o-visualization "$D/rarefaction_grouped.qzv"
    else
        log "No GROUP_COL or no metadata columns yet: drawing per-sample curves only."
    fi
    run qiime diversity alpha-rarefaction \
        --i-table "$D/table_taxfilt.qza" \
        --p-max-depth "$MEDIAN" --p-steps 20 \
        --p-metrics observed_features shannon \
        --o-visualization "$D/rarefaction_per_sample.qzv"
fi

need MIN_SAMPLE_DEPTH "Inspect $D/rarefaction_*.qzv and the Interactive Sample Detail tab of $D/table_taxfilt_summary.qzv; set the depth below which samples are considered failed (use 1 to keep all)."

have "$D/table_final.qza" || run qiime feature-table filter-samples \
    --i-table "$D/table_taxfilt.qza" \
    --p-min-frequency "$MIN_SAMPLE_DEPTH" \
    --o-filtered-table "$D/table_final.qza"

# Features that only occurred in removed samples are dropped from the sequences too.
have "$D/rep_seqs_final.qza" || run qiime feature-table filter-seqs \
    --i-data "$SEQS" \
    --i-table "$D/table_final.qza" \
    --o-filtered-data "$D/rep_seqs_final.qza"

summarize "$D/table_final.qza" "$D/table_final"

log "CHECKPOINT: in $D/table_final_summary.qzv choose SAMPLING_DEPTH (Interactive Sample"
log "Detail; group by ${GROUP_COL:-your main column} to see losses per group)."
