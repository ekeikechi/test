#!/usr/bin/env bash
# Stage 09: taxonomic composition and differential abundance.
#   Barplots    : taxa barplot (+ per-group sums)   [all sources] [MH_qiime2-amplicon-2024.5 §7]
#   Diff. abund.: ANCOM-BC2 on ASVs and at TAXA_COLLAPSE_LEVEL [gut_to_soil_tutorial]
#                 (classic ANCOM [MH_*] and ANCOM-BC [moving_pictures_tutorial] are older
#                 methods; ANCOM-BC2 adds bias correction, random effects and a
#                 pseudo-count sensitivity analysis)
# DA runs on counts that were NOT rarefied (table_final), as ANCOM-BC2 expects.
source "$(dirname "$0")/lib.sh"
D="${OUT}/09_composition_da"; mkdir -p "$D"
show_settings "$D" \
    Q2_EPOCH GROUP_COL DA_SAMPLE_FILTER DA_REFERENCE TAXONOMY_RANKS \
    TAXA_COLLAPSE_LEVEL NCORES
TABLE="${OUT}/06_filter_taxa_samples/table_final.qza"
TAX="${OUT}/05_taxonomy/taxonomy.qza"

have "$D/taxa_barplot.qzv" || run qiime taxa barplot \
    --i-table "$TABLE" \
    --i-taxonomy "$TAX" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$D/taxa_barplot.qzv"

if [[ -n "$GROUP_COL" ]]; then
    # Summed per group; useful for overview figures, not for statistics.
    have "$D/taxa_barplot_by_${GROUP_COL}.qzv" || {
        run qiime feature-table group \
            --i-table "$TABLE" \
            --p-axis sample \
            --p-mode sum \
            --m-metadata-file "$METADATA" \
            --m-metadata-column "$GROUP_COL" \
            --o-grouped-table "$D/table_by_${GROUP_COL}.qza"
        run qiime taxa barplot \
            --i-table "$D/table_by_${GROUP_COL}.qza" \
            --i-taxonomy "$TAX" \
            --o-visualization "$D/taxa_barplot_by_${GROUP_COL}.qzv"
    }
fi

need GROUP_COL "Set the metadata column to test for differential abundance."

DA_TABLE="$TABLE"
if [[ -n "$DA_SAMPLE_FILTER" ]]; then
    have "$D/table_da.qza" || run qiime feature-table filter-samples \
        --i-table "$TABLE" \
        --m-metadata-file "$METADATA" \
        --p-where "$DA_SAMPLE_FILTER" \
        --o-filtered-table "$D/table_da.qza"
    DA_TABLE="$D/table_da.qza"
fi

# Collapsed-table DA at the configured rank (number in TAXONOMY_RANKS; R uses the same).
RANK="$(rank_name "$TAXA_COLLAPSE_LEVEL" | tr '[:upper:]' '[:lower:]')"
log "Collapsing to taxonomy level ${TAXA_COLLAPSE_LEVEL} (${RANK}) of: ${TAXONOMY_RANKS}"
have "$D/table_da_${RANK}.qza" || run qiime taxa collapse \
    --i-table "$DA_TABLE" \
    --i-taxonomy "$TAX" \
    --p-level "$TAXA_COLLAPSE_LEVEL" \
    --o-collapsed-table "$D/table_da_${RANK}.qza"

ref_args=()
[[ -n "$DA_REFERENCE" ]] && ref_args=(--p-reference-levels "$DA_REFERENCE")

# Visualizer name changed: 2026.1 has ancombc2-visualizer; 2026.4+ uses da-barplot.
if epoch_ge 2026.4; then VIZ=da-barplot; else VIZ=ancombc2-visualizer; fi

for level in asv "$RANK"; do
    if [[ "$level" == "asv" ]]; then t="$DA_TABLE"; tax_args=(--i-taxonomy "$TAX"); else t="$D/table_da_${RANK}.qza"; tax_args=(); fi
    have "$D/ancombc2_${level}.qza" || run qiime composition ancombc2 \
        --i-table "$t" \
        --m-metadata-file "$METADATA" \
        --p-fixed-effects-formula "$GROUP_COL" \
        ${ref_args[@]+"${ref_args[@]}"} \
        --p-num-processes "$NCORES" \
        --o-ancombc2-output "$D/ancombc2_${level}.qza"
    have "$D/ancombc2_${level}.qzv" || run qiime composition "$VIZ" \
        --i-data "$D/ancombc2_${level}.qza" \
        ${tax_args[@]+"${tax_args[@]}"} \
        --o-visualization "$D/ancombc2_${level}.qzv"
    have "$D/ancombc2_${level}_table.qzv" || run qiime composition tabulate \
        --i-data "$D/ancombc2_${level}.qza" \
        --o-visualization "$D/ancombc2_${level}_table.qzv"
done

log "Treat significant features as hypotheses [gut_to_soil_tutorial]; confirm with the"
log "multi-method consensus in R (R/downstream_analysis.R) [MH_2026.3-R]."
