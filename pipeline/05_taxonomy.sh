#!/usr/bin/env bash
# Stage 05: taxonomic classification with a pre-trained Naive Bayes classifier.
# Sources: [MH_2026.3-workflow §3] [MH_qiime2-amplicon-2024.5 §3]
#          [gut_to_soil_tutorial: tabulate-seqs with taxonomy] [conceptual_overview: weighted classifiers]
# The classifier must have been trained with the scikit-learn version of this QIIME 2
# release, otherwise classify-sklearn stops with a version-mismatch ValueError.
source "$(dirname "$0")/lib.sh"
D="${OUT}/05_taxonomy"; mkdir -p "$D"
show_settings "$D" \
    CLASSIFIER CLASSIFIER_LABEL CLASSIFY_CONFIDENCE TAXONOMY_RANKS NCORES
SEQS="${OUT}/04_filter_features/rep_seqs.qza"

need CLASSIFIER "Download a classifier matching QIIME 2 ${Q2_EPOCH} (see docs/02-optimized-amplicon-sop.md, stage 5)."
if [[ ! -f "$CLASSIFIER" ]]; then
    [[ "$BUILD_PR2_CLASSIFIER" == "yes" ]] && die "Classifier not found: $CLASSIFIER. Run 05a_build_pr2_classifier.sh first."
    die "Classifier not found: $CLASSIFIER"
fi

# Reads are in amplicon orientation after primer trimming (and after denoise-ccs),
# so 'same' skips the orientation guess [MH_2026.3-workflow §3.2].
have "$D/taxonomy.qza" || run qiime feature-classifier classify-sklearn \
    --i-reads "$SEQS" \
    --i-classifier "$CLASSIFIER" \
    --p-confidence "$CLASSIFY_CONFIDENCE" \
    --p-read-orientation same \
    --p-n-jobs "$NCORES" \
    --o-classification "$D/taxonomy.qza"

# One table with sequence, taxonomy, frequency and BLAST links, for spot checks.
# --i-taxonomy takes "label:file"; the label must be a Python identifier, so
# "PR2_5.1.1" is passed as "PR2_5_1_1".
TAX_KEY="$(printf '%s' "${CLASSIFIER_LABEL:-taxonomy}" | sed -e 's/[^A-Za-z0-9_]/_/g' -e 's/^\([0-9]\)/_\1/')"
have "$D/rep_seqs_taxonomy.qzv" || run qiime feature-table tabulate-seqs \
    --i-data "$SEQS" \
    --i-taxonomy "${TAX_KEY:-taxonomy}:$D/taxonomy.qza" \
    --m-metadata-file "${OUT}/04_filter_features/table_feature_frequencies.qza" \
    --o-visualization "$D/rep_seqs_taxonomy.qzv"

log "CHECKPOINT: in $D/rep_seqs_taxonomy.qzv BLAST ~5 abundant ASVs from different phyla"
log "and compare with the classifier labels [MH_*] [moving_pictures_tutorial]."
log "Many ASVs unassigned at phylum level can mean wrong orientation, wrong classifier"
log "for the marker, or non-target amplification."
