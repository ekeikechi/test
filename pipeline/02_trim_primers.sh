#!/usr/bin/env bash
# Stage 02: remove primers with cutadapt (Illumina only).
# Sources: [MH_qiime2-amplicon-2024.5 §1.6-1.7] [MH_2026.3-workflow §2.3.1-2.3.2]
#          [experienced_researchers: removing non-biological sequences]
# PacBio primers are removed inside dada2 denoise-ccs (stage 03).
source "$(dirname "$0")/lib.sh"
D="${OUT}/02_trim"; mkdir -p "$D"
show_settings "$D" \
    PLATFORM LAYOUT MARKER PRIMER_F PRIMER_R NCORES
IN="${OUT}/01_import/reads.qza"

if [[ "$PLATFORM" == "pacbio" ]]; then
    log "PacBio: primers are handled by dada2 denoise-ccs; nothing to do."; exit 0
fi

if [[ -z "$PRIMER_F" && -z "$PRIMER_R" ]]; then
    log "No primers configured; linking the untrimmed reads forward."
    ln -sf "$IN" "$D/reads_trimmed.qza"
    exit 0
fi

# --p-discard-untrimmed drops reads without the primer (off-target / primer dimers);
# --p-no-indels matches MH practice for anchored primer sequences.
if ! have "$D/reads_trimmed.qza"; then
    if [[ "$LAYOUT" == "paired" ]]; then
        run qiime cutadapt trim-paired \
            --i-demultiplexed-sequences "$IN" \
            --p-front-f "$PRIMER_F" \
            --p-front-r "$PRIMER_R" \
            --p-discard-untrimmed \
            --p-no-indels \
            --p-cores "$NCORES" \
            --o-trimmed-sequences "$D/reads_trimmed.qza"
    else
        run qiime cutadapt trim-single \
            --i-demultiplexed-sequences "$IN" \
            --p-front "$PRIMER_F" \
            --p-discard-untrimmed \
            --p-no-indels \
            --p-cores "$NCORES" \
            --o-trimmed-sequences "$D/reads_trimmed.qza"
    fi
fi

have "$D/reads_trimmed_summary.qzv" || run qiime demux summarize \
    --i-data "$D/reads_trimmed.qza" \
    --o-visualization "$D/reads_trimmed_summary.qzv"

log "CHECKPOINT: compare $D/reads_trimmed_summary.qzv with 01_import/reads_summary.qzv."
log "Most reads should survive; a large loss usually means wrong primer sequences."
log "Then set DADA2_TRUNC_F/R (or DEBLUR_TRIM_LENGTH) from the quality plot."
