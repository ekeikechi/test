#!/usr/bin/env bash
# Stage 07: phylogeny for Faith PD / UniFrac.
#   16S : SEPP insertion into a reference tree, then drop unplaced ASVs
#         [MH_2026.3-workflow §5.1] [experienced_researchers: reference trees are higher quality]
#         (filter-features is the q2-fragment-insertion step for unplaced fragments)
#   18S : de novo MAFFT + FastTree + midpoint root, as one pipeline
#         [moving_pictures_tutorial] (same steps as MH's four separate commands)
#   ITS / other : no tree; use non-phylogenetic or kmer-based metrics [gut_to_soil_tutorial]
source "$(dirname "$0")/lib.sh"
D="${OUT}/07_phylogeny"; mkdir -p "$D"
show_settings "$D" \
    MARKER TREE_METHOD SEPP_REF NCORES
TABLE="${OUT}/06_filter_taxa_samples/table_final.qza"
SEQS="${OUT}/06_filter_taxa_samples/rep_seqs_final.qza"

METHOD="$TREE_METHOD"
if [[ "$METHOD" == "auto" ]]; then
    case "$MARKER" in
        16S) METHOD=sepp ;;
        18S) METHOD=denovo ;;
        *)   METHOD=none ;;
    esac
fi
log "Tree method: ${METHOD}"

case "$METHOD" in
    sepp)
        need SEPP_REF "Download a SEPP reference (e.g. Greengenes2 2022.10 backbone, or sepp-refs-silva-128.qza)."
        have "$D/tree.qza" || run qiime fragment-insertion sepp \
            --i-representative-sequences "$SEQS" \
            --i-reference-database "$SEPP_REF" \
            --p-threads "$NCORES" \
            --o-tree "$D/tree.qza" \
            --o-placements "$D/placements.qza"
        have "$D/table_placed.qza" || run qiime fragment-insertion filter-features \
            --i-table "$TABLE" \
            --i-tree "$D/tree.qza" \
            --o-filtered-table "$D/table_placed.qza" \
            --o-removed-table "$D/table_unplaced.qza"
        # The removed table is empty when every ASV was placed, and summarizing it can fail.
        summarize "$D/table_unplaced.qza" "$D/table_unplaced" \
            || log "No summary for unplaced ASVs (likely none were rejected)."
        log "CHECKPOINT: $D/table_unplaced_summary.qzv lists reads lost because SEPP could not"
        log "place their ASVs (usually few). Phylogenetic metrics use table_placed.qza." ;;
    denovo)
        have "$D/tree.qza" || run qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences "$SEQS" \
            --p-n-threads "$NCORES" \
            --o-alignment "$D/aligned_rep_seqs.qza" \
            --o-masked-alignment "$D/masked_aligned_rep_seqs.qza" \
            --o-tree "$D/unrooted_tree.qza" \
            --o-rooted-tree "$D/tree.qza"
        ln -sf "$TABLE" "$D/table_placed.qza" ;;
    none)
        log "No tree for MARKER=${MARKER}; stage 08 will use non-phylogenetic metrics." ;;
    *) die "TREE_METHOD must be auto, sepp, denovo or none" ;;
esac
