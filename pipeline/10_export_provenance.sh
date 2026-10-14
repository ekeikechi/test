#!/usr/bin/env bash
# Stage 10: export plain-text files for R / other tools, and record provenance.
#   Export     : FASTA, BIOM (+ taxonomy), TSV, Newick  [MH_qiime2-amplicon-2024.5 §10] [MH_2026.3-workflow §6]
#   Provenance : replay-provenance / replay-supplement   [gut_to_soil_tutorial: provenance replay]
source "$(dirname "$0")/lib.sh"
D="${OUT}/10_export"; mkdir -p "$D"
show_settings "$D" \
    TAXONOMY_RANKS TAXA_COLLAPSE_LEVEL OVERVIEW_RANK_LEVEL SAMPLING_DEPTH GROUP_COL \
    DA_REFERENCE RANDOM_SEED SAMPLE_HABITAT TRAIT_MIN_AGREEMENT BUILD_PR2_CLASSIFIER PR2_VERSION \
    MIN_SAMPLE_DEPTH DEPTH_SENSITIVITY_GROUPS DEPTH_SENSITIVITY_MIN_KEEP
TABLE="${OUT}/06_filter_taxa_samples/table_final.qza"
SEQS="${OUT}/06_filter_taxa_samples/rep_seqs_final.qza"
TAX="${OUT}/05_taxonomy/taxonomy.qza"
TREE="${OUT}/07_phylogeny/tree.qza"

X="$D/qiime_export"   # raw `qiime tools export` folders, kept for reference
mkdir -p "$X"
export_once() { have "$X/$2" || run qiime tools export --input-path "$1" --output-path "$X/$2"; }

# Unrarefied counts: R does its own normalisation per analysis.
export_once "$TABLE" table
export_once "$TAX"   taxonomy
export_once "$SEQS"  rep_seqs
cp "$X/table/feature-table.biom"      "$D/feature-table.biom"
cp "$X/taxonomy/taxonomy.tsv"         "$D/taxonomy.tsv"
cp "$X/rep_seqs/dna-sequences.fasta"  "$D/rep_seqs.fasta"
cp "$METADATA"                        "$D/metadata.tsv"

have "$D/feature-table.tsv" \
    || run biom convert -i "$D/feature-table.biom" -o "$D/feature-table.tsv" --to-tsv

# BIOM/TSV with taxonomy as observation metadata (MH format). The header is renamed in
# a separate file so the exported taxonomy.tsv stays untouched.
sed -e '1 s/Feature ID/#Feature ID/' -e '1 s/Taxon/taxonomy/' \
    "$D/taxonomy.tsv" > "$X/taxonomy_biom_header.tsv"
have "$D/feature-table_w_tax.biom" || run biom add-metadata \
    -i "$D/feature-table.biom" \
    -o "$D/feature-table_w_tax.biom" \
    --observation-metadata-fp "$X/taxonomy_biom_header.tsv" \
    --sc-separated taxonomy
have "$D/feature-table_w_tax.txt" || run biom convert \
    -i "$D/feature-table_w_tax.biom" \
    -o "$D/feature-table_w_tax.txt" \
    --to-tsv --header-key taxonomy

if [[ -e "$TREE" ]]; then
    export_once "$TREE" tree
    cp "$X/tree/tree.nwk" "$D/tree.nwk"
fi

# The taxonomy must not be deeper than TAXONOMY_RANKS, or ranks would be mislabelled.
DEPTH="$(awk -F'\t' 'NR > 1 { sub(/;[[:space:]]*$/, "", $2); n = split($2, a, ";"); if (n > m) m = n }
                     END { print m + 0 }' "$D/taxonomy.tsv")"
(( DEPTH <= $(rank_count) )) \
    || die "taxonomy.tsv has up to ${DEPTH} ranks but TAXONOMY_RANKS lists $(rank_count): ${TAXONOMY_RANKS}"
(( DEPTH == $(rank_count) )) \
    || log "Note: the deepest assignment has ${DEPTH} of $(rank_count) ranks (fine if no ASV reached the last rank)."

# Single hand-off of shared settings to R (R/downstream_analysis.R reads this file),
# so QIIME 2 and R analyse the same rank, depth, group, reference and seed.
{
    printf 'setting\tvalue\n'
    printf 'preset\t%s\n'           "$PRESET"
    printf 'taxonomy_ranks\t%s\n'   "$TAXONOMY_RANKS"
    printf 'collapse_level\t%s\n'   "$TAXA_COLLAPSE_LEVEL"
    printf 'overview_level\t%s\n'   "$OVERVIEW_RANK_LEVEL"
    printf 'sampling_depth\t%s\n'   "$SAMPLING_DEPTH"
    printf 'min_sample_depth\t%s\n' "$MIN_SAMPLE_DEPTH"
    printf 'depth_sensitivity_groups\t%s\n'   "$DEPTH_SENSITIVITY_GROUPS"
    printf 'depth_sensitivity_min_keep\t%s\n' "$DEPTH_SENSITIVITY_MIN_KEEP"
    printf 'group_col\t%s\n'        "$GROUP_COL"
    printf 'da_reference\t%s\n'     "$DA_REFERENCE"
    printf 'random_seed\t%s\n'      "$RANDOM_SEED"
    # read by R/pr2_traits.R
    printf 'sample_habitat\t%s\n'      "$SAMPLE_HABITAT"
    printf 'trait_min_agreement\t%s\n' "$TRAIT_MIN_AGREEMENT"
    printf 'pr2_version\t%s\n'         "$([[ "$BUILD_PR2_CLASSIFIER" == yes ]] && echo "$PR2_VERSION")"
} > "$D/analysis_settings.tsv"

log "Provenance: a re-runnable script and a methods supplement (citations + script)"
have "$D/replayed_commands.sh" || run qiime tools replay-provenance \
    --in-fp "$OUT" --recurse \
    --out-fp "$D/replayed_commands.sh"
have "$D/reproducibility_supplement.zip" || run qiime tools replay-supplement \
    --in-fp "$OUT" --recurse \
    --out-fp "$D/reproducibility_supplement.zip"

log "Done. R inputs are in $D (feature-table.tsv, taxonomy.tsv, metadata.tsv, tree.nwk)."
if [[ "$BUILD_PR2_CLASSIFIER" == yes ]]; then
    log "PR2 traits next: Rscript R/pr2_traits.R $D ${OUT}/11_traits"
fi
