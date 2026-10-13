#!/usr/bin/env bash
# Stage 05a (optional): build a PR2 Naive Bayes classifier with RESCRIPt.
#   get-pr2-data -> dereplicate -> fit-classifier-naive-bayes
# Training in the same QIIME 2 environment guarantees the scikit-learn version matches
# classify-sklearn in stage 05. The full-length reference is used, as the SOP recommends
# for the other databases [MH_2026.3-workflow] [gut_to_soil_tutorial: RESCRIPt].
# Runs only when BUILD_PR2_CLASSIFIER=yes (set by PRESET=18S_PR2).
source "$(dirname "$0")/lib.sh"
D="${OUT}/05a_pr2_classifier"; mkdir -p "$D"
show_settings "$D" \
    BUILD_PR2_CLASSIFIER PR2_VERSION CLASSIFIER TAXONOMY_RANKS NCORES

if [[ "$BUILD_PR2_CLASSIFIER" != "yes" ]]; then
    log "BUILD_PR2_CLASSIFIER is not 'yes'; nothing to do."; exit 0
fi
TARGET="${D}/pr2_${PR2_VERSION}_classifier.qza"
if [[ "$CLASSIFIER" != "$TARGET" ]]; then
    log "Note: stage 05 will use CLASSIFIER=${CLASSIFIER:-<empty>}, not the classifier built here (${TARGET})."
fi
# The PR2 taxonomy RESCRIPt produces always has these 9 ranks.
[[ "$(rank_count)" -eq 9 ]] \
    || die "PR2 has 9 ranks but TAXONOMY_RANKS lists $(rank_count): ${TAXONOMY_RANKS}"

# Labels are prefixed d__ sgr__ dv__ dvs__ c__ o__ f__ g__ s__. PR2 also contains
# Bacteria/Archaea outgroups and organelle sequences (:plas, :mito, ...); they stay in
# the reference so such reads are recognised, then removed in stage 06 (TAXA_EXCLUDE).
have "$D/pr2_sequences.qza" || run qiime rescript get-pr2-data \
    --p-version "$PR2_VERSION" \
    --o-pr2-sequences "$D/pr2_sequences.qza" \
    --o-pr2-taxonomy "$D/pr2_taxonomy.qza"

# Collapse identical reference sequences. --p-rank-handles disable stops RESCRIPt from
# backfilling missing ranks with its default 7 NCBI-style prefixes, which would not
# match PR2's 9 ranks.
have "$D/pr2_derep_sequences.qza" || run qiime rescript dereplicate \
    --i-sequences "$D/pr2_sequences.qza" \
    --i-taxa "$D/pr2_taxonomy.qza" \
    --p-mode uniq \
    --p-rank-handles disable \
    --p-threads "$NCORES" \
    --o-dereplicated-sequences "$D/pr2_derep_sequences.qza" \
    --o-dereplicated-taxa "$D/pr2_derep_taxonomy.qza"

# Memory-hungry: expect it to need a large-memory node for the full PR2 reference
# (not measured here); see the SOP for lowering memory use.
have "$TARGET" || run qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$D/pr2_derep_sequences.qza" \
    --i-reference-taxonomy "$D/pr2_derep_taxonomy.qza" \
    --o-classifier "$TARGET"

log "Classifier ready: ${TARGET}"
