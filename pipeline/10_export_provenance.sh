#!/usr/bin/env bash
# Stage 10: export plain-text files for R / other tools, and record provenance.
#   Export     : FASTA, BIOM (+ taxonomy), TSV, Newick  [MH_qiime2-amplicon-2024.5 §10] [MH_2026.3-workflow §6]
#   Provenance : replay-provenance / replay-supplement   [gut_to_soil_tutorial: provenance replay]
source "$(dirname "$0")/lib.sh"
D="${OUT}/10_export"; mkdir -p "$D"
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

log "Provenance: a re-runnable script and a methods supplement (citations + script)"
have "$D/replayed_commands.sh" || run qiime tools replay-provenance \
    --in-fp "$OUT" --recurse \
    --out-fp "$D/replayed_commands.sh"
have "$D/reproducibility_supplement.zip" || run qiime tools replay-supplement \
    --in-fp "$OUT" --recurse \
    --out-fp "$D/reproducibility_supplement.zip"

log "Done. R inputs are in $D (feature-table.tsv, taxonomy.tsv, metadata.tsv, tree.nwk)."
