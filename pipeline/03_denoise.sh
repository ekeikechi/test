#!/usr/bin/env bash
# Stage 03: denoise reads into ASVs.
#   Illumina default : DADA2 denoise-paired / denoise-single   [gut_to_soil_tutorial] [moving_pictures_tutorial]
#   Illumina option  : vsearch merge -> q-score -> Deblur      [MH_qiime2-amplicon-2024.5 §2] [MH_2026.3-workflow §2.3.3]
#   PacBio HiFi      : DADA2 denoise-ccs                       [MH_2026.3-workflow §2.3 PacBio]
# All branches write D/table.qza and D/rep_seqs.qza so later stages need no branching.
source "$(dirname "$0")/lib.sh"
D="${OUT}/03_denoise"; mkdir -p "$D"
show_settings "$D" \
    PLATFORM LAYOUT MARKER DENOISER DADA2_TRIM_LEFT_F DADA2_TRIM_LEFT_R \
    DADA2_TRUNC_F DADA2_TRUNC_R DADA2_TRUNC_MIN_MEDIAN_Q DADA2_TRUNC_MIN_READ_FRACTION \
    DADA2_TRUNC_SUM_TARGET DADA2_MAX_EE_F DADA2_MAX_EE_R DADA2_POOLING \
    DEBLUR_TRIM_LENGTH DEBLUR_REF_SEQS \
    CCS_FRONT CCS_ADAPTER CCS_MIN_LEN CCS_MAX_LEN CCS_MAX_EE NCORES

if [[ "$PLATFORM" == "pacbio" ]]; then
    IN="${OUT}/01_import/reads.qza"
else
    IN="${OUT}/02_trim/reads_trimmed.qza"
fi
[[ -e "$IN" ]] || die "Input not found: $IN (run the previous stage first)."

# auto_truncation SUMMARY_QZV: replace DADA2_TRUNC_F/R ("auto") with lengths chosen by
# auto_truncation.py from the trimmed reads' per-position quality; see config.sh.
auto_truncation() {
    local qzv="$1" sumdir="${OUT}/02_trim/reads_trimmed_summary" dec="$D/truncation_auto.tsv"
    local picked args
    [[ "$LAYOUT" != "paired" || ( "$DADA2_TRUNC_F" == auto && "$DADA2_TRUNC_R" == auto ) ]] \
        || die "Set both DADA2_TRUNC_F and DADA2_TRUNC_R to \"auto\", or both to numbers."
    [[ -e "$qzv" ]] || die "Automatic truncation needs ${qzv} (run stage 02 first)."
    have "$sumdir/forward-seven-number-summaries.tsv" || run qiime tools export \
        --input-path "$qzv" --output-path "$sumdir"
    args=(--forward "$sumdir/forward-seven-number-summaries.tsv"
          --min-median-q "$DADA2_TRUNC_MIN_MEDIAN_Q"
          --min-read-fraction "$DADA2_TRUNC_MIN_READ_FRACTION" --out "$dec")
    if [[ "$LAYOUT" == "paired" ]]; then
        args+=(--reverse "$sumdir/reverse-seven-number-summaries.tsv")
        if [[ -n "$DADA2_TRUNC_SUM_TARGET" ]]; then
            args+=(--sum-target "$DADA2_TRUNC_SUM_TARGET")
        else
            log "WARNING: DADA2_TRUNC_SUM_TARGET is empty, so the automatic truncation does not"
            log "check that read pairs can still merge. Check the merged column of denoising_stats."
        fi
    fi
    log "Choosing truncation lengths automatically (quality rule; decisions in ${dec})"
    picked="$(python3 "${PIPELINE_DIR}/auto_truncation.py" "${args[@]}")" \
        || die "Automatic truncation stopped; see the message above and ${dec}."
    read -r DADA2_TRUNC_F DADA2_TRUNC_R <<< "$picked"
    awk -F'\t' 'NR > 1 { printf "  %-32s %s\n", $1, $2 }' "$dec" >&2
    # record the chosen lengths in place of "auto" in this stage's settings report
    awk -F'\t' -v OFS='\t' -v f="$DADA2_TRUNC_F" -v r="${DADA2_TRUNC_R:-}" '
        ($1 == "DADA2_TRUNC_F" || ($1 == "DADA2_TRUNC_R" && r != "")) && $2 == "auto" {
            $3 = "auto, chosen in truncation_auto.tsv; " $3; $2 = ($1 == "DADA2_TRUNC_F" ? f : r) }
        { print }' "$D/settings_used.tsv" > "$D/settings_used.tsv.new"
    mv "$D/settings_used.tsv.new" "$D/settings_used.tsv"
    log "Automatic truncation: DADA2_TRUNC_F=${DADA2_TRUNC_F}${DADA2_TRUNC_R:+ DADA2_TRUNC_R=${DADA2_TRUNC_R}}"
}

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
    need DADA2_TRUNC_F "Choose it from the Interactive Quality Plot in ${TRIMMED_QZV}, or use \"auto\"."
    [[ "$LAYOUT" == "paired" ]] && need DADA2_TRUNC_R "Choose it from the Interactive Quality Plot in ${TRIMMED_QZV}, or use \"auto\"."
    if [[ "$DADA2_TRUNC_F" == auto || ( "$LAYOUT" == "paired" && "$DADA2_TRUNC_R" == auto ) ]]; then
        auto_truncation "$TRIMMED_QZV"
    fi
    if [[ "$LAYOUT" == "paired" ]]; then
        log "DADA2 denoise-paired (merging, chimera removal and phiX filtering are built in)"
        log "Reminder: trunc-len-f + trunc-len-r must exceed the trimmed amplicon length by >= 12 nt."
        have "$D/table.qza" || run qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$IN" \
            --p-trim-left-f "$DADA2_TRIM_LEFT_F" --p-trim-left-r "$DADA2_TRIM_LEFT_R" \
            --p-trunc-len-f "$DADA2_TRUNC_F" --p-trunc-len-r "$DADA2_TRUNC_R" \
            --p-max-ee-f "$DADA2_MAX_EE_F" --p-max-ee-r "$DADA2_MAX_EE_R" \
            --p-pooling-method "$DADA2_POOLING" \
            --p-n-threads "$NCORES" \
            "${dada2_outputs[@]}"
    else
        log "DADA2 denoise-single"
        have "$D/table.qza" || run qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$IN" \
            --p-trim-left "$DADA2_TRIM_LEFT_F" \
            --p-trunc-len "$DADA2_TRUNC_F" \
            --p-max-ee "$DADA2_MAX_EE_F" \
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
