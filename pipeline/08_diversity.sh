#!/usr/bin/env bash
# Stage 08: alpha/beta diversity with repeated rarefaction, then group statistics.
#   Engine : q2-boots core-metrics, rarefaction (no replacement), n iterations, averaged
#            [gut_to_soil_tutorial] [conceptual_overview: q2-boots]; replaces the single
#            rarefaction of core-metrics-phylogenetic [MH_*] [moving_pictures_tutorial].
#   Tests  : alpha-group-significance, alpha-correlation, beta-group-significance (pairwise)
#            [moving_pictures_tutorial] [MH_2026.3-stats-visualization]; PERMDISP added
#            because PERMANOVA is sensitive to unequal dispersion; adonis for multi-factor.
source "$(dirname "$0")/lib.sh"
D="${OUT}/08_diversity"; mkdir -p "$D"
show_settings "$D" \
    Q2_EPOCH SAMPLING_DEPTH DEPTH_SENSITIVITY_GROUPS DEPTH_SENSITIVITY_MIN_KEEP \
    BOOTS_N RANDOM_SEED RUN_KMER_DIVERSITY GROUP_COL \
    ADONIS_FORMULA NCORES
TREE="${OUT}/07_phylogeny/tree.qza"

need SAMPLING_DEPTH "Choose it in 06_filter_taxa_samples/table_final_summary.qzv (Interactive Sample Detail), or use \"auto\"."

if [[ -e "$TREE" ]]; then
    TABLE="${OUT}/07_phylogeny/table_placed.qza"
    phylo_args=(--i-phylogeny "$TREE")
else
    TABLE="${OUT}/06_filter_taxa_samples/table_final.qza"
    phylo_args=()
fi
seed_args=()
epoch_ge 2026.4 && seed_args=(--p-random-seed "$RANDOM_SEED")

# Depth rules (pipeline/choose_depth.py): "auto" main depth and the sensitivity depth.
SENS_DEPTH=""
if [[ "$SAMPLING_DEPTH" == auto || -n "$DEPTH_SENSITIVITY_GROUPS" ]]; then
    sample_depths "$TABLE" "$D/table_biom" "$D/sample_depths.tsv"
    grp_args=()
    [[ -n "$DEPTH_SENSITIVITY_GROUPS" ]] && grp_args=(--metadata "$METADATA"
        --groups "$DEPTH_SENSITIVITY_GROUPS" --min-keep "$DEPTH_SENSITIVITY_MIN_KEEP"
        --groups-out "$D/depth_by_group.tsv")
    picked="$(python3 "${PIPELINE_DIR}/choose_depth.py" --depths "$D/sample_depths.tsv" \
        --out "$D/depth_choice.tsv" ${grp_args[@]+"${grp_args[@]}"})" \
        || die "Depth choice failed; see the message above."
    read -r MAIN_DEPTH SENS_DEPTH <<< "$picked"
    awk -F'\t' 'NR > 1 { printf "  %-30s %s\n", $1, $2 }' "$D/depth_choice.tsv" >&2
    if [[ "$SAMPLING_DEPTH" == auto ]]; then
        SAMPLING_DEPTH="$MAIN_DEPTH"
        awk -F'\t' -v OFS='\t' -v d="$SAMPLING_DEPTH" \
            '$1 == "SAMPLING_DEPTH" { $3 = "auto, chosen in depth_choice.tsv; " $3; $2 = d } { print }' \
            "$D/settings_used.tsv" > "$D/settings_used.tsv.new" && mv "$D/settings_used.tsv.new" "$D/settings_used.tsv"
    fi
    log "Rarefaction depth ${SAMPLING_DEPTH}${SENS_DEPTH:+; sensitivity depth ${SENS_DEPTH}}"
fi

# core_metrics OUTDIR DEPTH: q2-boots core-metrics; collection outputs are written as
# directories: <dir>/<metric>.qza (+ .order).
core_metrics() {
    local dir="$1" depth="$2"
    have "$dir/alpha" && return 0
    mkdir -p "$dir"
    run qiime boots core-metrics \
        --i-table "$TABLE" \
        ${phylo_args[@]+"${phylo_args[@]}"} \
        --m-metadata-file "$METADATA" \
        --p-sampling-depth "$depth" \
        --p-n "$BOOTS_N" \
        --p-no-replacement \
        --p-alpha-average-method median \
        --p-beta-average-method medoid \
        ${seed_args[@]+"${seed_args[@]}"} \
        --o-resampled-tables "$dir/resampled_tables" \
        --o-alpha-diversities "$dir/alpha" \
        --o-distance-matrices "$dir/distance" \
        --o-pcoas "$dir/pcoa" \
        --o-emperor-plots "$dir/emperor" \
        --o-scatter-plot "$dir/scatter_plot.qzv"
}
core_metrics "$D/core" "$SAMPLING_DEPTH"
if [[ -n "$SENS_DEPTH" && "$SENS_DEPTH" != "$SAMPLING_DEPTH" ]]; then
    log "Sensitivity run at ${SENS_DEPTH} reads (samples dropped per group: $D/depth_by_group.tsv)"
    core_metrics "$D/sensitivity/core" "$SENS_DEPTH"
fi

log "Alpha diversity tests (Kruskal-Wallis for categories, Spearman for numeric columns)"
for vec in "$D"/core/alpha/*.qza; do
    m="$(basename "$vec" .qza)"
    have "$D/alpha_${m}_groups.qzv" || run qiime diversity alpha-group-significance \
        --i-alpha-diversity "$vec" \
        --m-metadata-file "$METADATA" \
        --o-visualization "$D/alpha_${m}_groups.qzv" \
        || log "alpha-group-significance skipped for ${m} (no categorical metadata columns?)"
    have "$D/alpha_${m}_correlation.qzv" || run qiime diversity alpha-correlation \
        --i-alpha-diversity "$vec" \
        --m-metadata-file "$METADATA" \
        --p-method spearman \
        --o-visualization "$D/alpha_${m}_correlation.qzv" \
        || log "alpha-correlation skipped for ${m} (no numeric metadata columns?)"
done

if [[ -n "$GROUP_COL" ]]; then
    log "Beta diversity: PERMANOVA (pairwise) and PERMDISP by ${GROUP_COL}"
    for dm in "$D"/core/distance/*.qza; do
        m="$(basename "$dm" .qza)"
        for method in permanova permdisp; do
            have "$D/beta_${m}_${GROUP_COL}_${method}.qzv" || run qiime diversity beta-group-significance \
                --i-distance-matrix "$dm" \
                --m-metadata-file "$METADATA" \
                --m-metadata-column "$GROUP_COL" \
                --p-method "$method" \
                --p-pairwise \
                --p-permutations 999 \
                --o-visualization "$D/beta_${m}_${GROUP_COL}_${method}.qzv"
        done
    done
else
    log "GROUP_COL not set; skipping beta-group-significance."
fi

if [[ -n "$ADONIS_FORMULA" ]]; then
    log "Multi-factor PERMANOVA (adonis2, sequential terms): ${ADONIS_FORMULA}"
    for dm in "$D"/core/distance/*.qza; do
        m="$(basename "$dm" .qza)"
        have "$D/adonis_${m}.qzv" || run qiime diversity adonis \
            --i-distance-matrix "$dm" \
            --m-metadata-file "$METADATA" \
            --p-formula "$ADONIS_FORMULA" \
            --p-permutations 999 \
            --p-n-jobs "$NCORES" \
            --o-visualization "$D/adonis_${m}.qzv"
    done
fi

if [[ "$RUN_KMER_DIVERSITY" == "yes" ]] && ! have "$D/kmer/scatter_plot.qzv"; then
    log "Tree-free kmer-based diversity [gut_to_soil_tutorial]"
    mkdir -p "$D/kmer"
    run qiime boots kmer-diversity \
        --i-table "${OUT}/06_filter_taxa_samples/table_final.qza" \
        --i-sequences "${OUT}/06_filter_taxa_samples/rep_seqs_final.qza" \
        --m-metadata-file "$METADATA" \
        --p-sampling-depth "$SAMPLING_DEPTH" \
        --p-n "$BOOTS_N" \
        --p-no-replacement \
        ${seed_args[@]+"${seed_args[@]}"} \
        --o-resampled-tables "$D/kmer/resampled_tables" \
        --o-kmer-tables "$D/kmer/kmer_tables" \
        --o-alpha-diversities "$D/kmer/alpha" \
        --o-distance-matrices "$D/kmer/distance" \
        --o-pcoas "$D/kmer/pcoa" \
        --o-scatter-plot "$D/kmer/scatter_plot.qzv"
fi

log "Open $D/core/scatter_plot.qzv and $D/core/emperor/*.qzv for ordinations."
log "Interpret PERMANOVA only alongside PERMDISP: a significant PERMDISP means group"
log "differences may reflect spread rather than location."
