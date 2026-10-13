#!/usr/bin/env bash
# Stage 03: denoise reads into ASVs.
#   Illumina default : DADA2 denoise-paired / denoise-single   [gut_to_soil_tutorial] [moving_pictures_tutorial]
#   Illumina option  : vsearch merge -> q-score -> Deblur      [MH_qiime2-amplicon-2024.5 §2] [MH_2026.3-workflow §2.3.3]
#   PacBio HiFi      : DADA2 denoise-ccs                       [MH_2026.3-workflow §2.3 PacBio]
# All branches write D/table.qza and D/rep_seqs.qza so later stages need no branching.
source "$(dirname "$0")/lib.sh"
D="${OUT}/03_denoise"; mkdir -p "$D"

if [[ "$PLATFORM" == "pacbio" ]]; then
    IN="${OUT}/01_import/reads.qza"
else
    IN="${OUT}/02_trim/reads_trimmed.qza"
fi
[[ -e "$IN" ]] || die "Input not found: $IN (run the previous stage first)."

dada2_outputs=(
    --o-table "$D/table.qza"
    --o-representative-sequences "$D/rep_seqs.qza"
    --o-denoising-stats "$D/denoising_stats.qza"
    --o-base-transition-stats "$D/base_transition_stats.qza"
)

if [[ "$PLATFORM" == "pacbio" ]]; then
    log "DADA2 denoise-ccs (primer removal, orientation and length gating inside DADA2)"
    have "$D/table.qza" || run qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs "$IN" \
        --p-front "$CCS_FRONT" \
        --p-adapter "$CCS_ADAPTER" \
        --p-min-len "$CCS_MIN_LEN" --p-max-len "$CCS_MAX_LEN" \
        --p-max-ee "$CCS_MAX_EE" \
        --p-pooling-method pseudo \
        --p-n-reads-learn 100000 \
        --p-n-threads "$NCORES" \
        "${dada2_outputs[@]}" \
        --verbose

elif [[ "$DENOISER" == "dada2" ]]; then
    TRIMMED_QZV="${OUT}/02_trim/reads_trimmed_summary.qzv"
    need DADA2_TRUNC_F "Choose it from the Interactive Quality Plot in ${TRIMMED_QZV}."
    if [[ "$LAYOUT" == "paired" ]]; then
        need DADA2_TRUNC_R "Choose it from the Interactive Quality Plot in ${TRIMMED_QZV}."
        log "DADA2 denoise-paired (merging, chimera removal and phiX filtering are built in)"
        log "Reminder: trunc-len-f + trunc-len-r must exceed the trimmed amplicon length by >= 12 nt."
        have "$D/table.qza" || run qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$IN" \
            --p-trim-left-f "$DADA2_TRIM_LEFT_F" --p-trim-left-r "$DADA2_TRIM_LEFT_R" \
            --p-trunc-len-f "$DADA2_TRUNC_F" --p-trunc-len-r "$DADA2_TRUNC_R" \
            --p-pooling-method "$DADA2_POOLING" \
            --p-n-threads "$NCORES" \
            "${dada2_outputs[@]}"
    else
        log "DADA2 denoise-single"
        have "$D/table.qza" || run qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$IN" \
            --p-trim-left "$DADA2_TRIM_LEFT_F" \
            --p-trunc-len "$DADA2_TRUNC_F" \
            --p-pooling-method "$DADA2_POOLING" \
            --p-n-threads "$NCORES" \
            "${dada2_outputs[@]}"
    fi

elif [[ "$DENOISER" == "deblur" ]]; then
    QF_IN="$IN"
    if [[ "$LAYOUT" == "paired" ]]; then
        log "Merging read pairs with vsearch (Deblur needs single reads)"
        have "$D/reads_merged.qza" || run qiime vsearch merge-pairs \
            --i-demultiplexed-seqs "$IN" \
            --p-threads "$NCORES" \
            --o-merged-sequences "$D/reads_merged.qza" \
            --o-unmerged-sequences "$D/reads_unmerged.qza"
        QF_IN="$D/reads_merged.qza"
    fi
    log "Quality filtering (recommended before Deblur)"
    have "$D/reads_qfilt.qza" || run qiime quality-filter q-score \
        --i-demux "$QF_IN" \
        --o-filtered-sequences "$D/reads_qfilt.qza" \
        --o-filter-stats "$D/qfilt_stats.qza"
    have "$D/reads_qfilt_summary.qzv" || run qiime demux summarize \
        --i-data "$D/reads_qfilt.qza" \
        --o-visualization "$D/reads_qfilt_summary.qzv"

    need DEBLUR_TRIM_LENGTH "Read the 'sequence length summary' percentiles in $D/reads_qfilt_summary.qzv (check that merging/filtering kept most reads)."
    # --p-min-reads 1 keeps rare ASVs so that rare-feature filtering is applied the
    # same way for every denoiser in stage 04 (Deblur's default would be 10).
    deblur_args=(
        --i-demultiplexed-seqs "$D/reads_qfilt.qza"
        --p-trim-length "$DEBLUR_TRIM_LENGTH"
        --p-sample-stats
        --p-min-reads 1
        --p-jobs-to-start "$NCORES"
        --o-table "$D/table.qza"
        --o-representative-sequences "$D/rep_seqs.qza"
        --o-stats "$D/deblur_stats.qza"
    )
    if [[ "$MARKER" == "16S" ]]; then
        have "$D/table.qza" || run qiime deblur denoise-16S "${deblur_args[@]}"
    else
        need DEBLUR_REF_SEQS "Deblur on ${MARKER} needs a positive-filter reference (e.g. PR2 for 18S, UNITE for ITS) imported as FeatureData[Sequence]."
        have "$D/table.qza" || run qiime deblur denoise-other \
            --i-reference-seqs "$DEBLUR_REF_SEQS" "${deblur_args[@]}"
    fi
    have "$D/deblur_stats.qzv" || run qiime deblur visualize-stats \
        --i-deblur-stats "$D/deblur_stats.qza" \
        --o-visualization "$D/deblur_stats.qzv"
else
    die "DENOISER must be dada2 or deblur"
fi

if [[ -e "$D/denoising_stats.qza" ]] && ! have "$D/denoising_stats.qzv"; then
    # Viewing the stats next to the sample metadata shows which sample types lose reads.
    run qiime metadata tabulate \
        --m-input-file "$METADATA" \
        --m-input-file "$D/denoising_stats.qza" \
        --o-visualization "$D/denoising_stats.qzv"
fi

summarize "$D/table.qza" "$D/table"
have "$D/rep_seqs.qzv" || run qiime feature-table tabulate-seqs \
    --i-data "$D/rep_seqs.qza" \
    --m-metadata-file "$D/table_feature_frequencies.qza" \
    --o-visualization "$D/rep_seqs.qzv"

log "CHECKPOINT: in $D/denoising_stats.qzv (or deblur_stats.qzv) check the fraction of"
log "input reads retained per sample. Large losses at filtering mean truncation is too"
log "aggressive; at merging, the overlap is too short; at chimera removal, primers remain."
