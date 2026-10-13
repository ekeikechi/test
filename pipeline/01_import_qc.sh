#!/usr/bin/env bash
# Stage 01: validate metadata, raw-read QC, import, and summarize.
# Sources: [MH_qiime2-amplicon-2024.5 §1] [MH_2026.3-workflow §1-2.2]
#          [gut_to_soil_tutorial: metadata tabulate, demux summarize]
#          [experienced_researchers: manifest import "will always work"]
source "$(dirname "$0")/lib.sh"
D="${OUT}/01_import"; mkdir -p "$D"
show_settings "$D" \
    Q2_EPOCH PRESET PLATFORM LAYOUT IMPORT_MODE RAW_DIR MANIFEST METADATA RUN_FASTQC \
    NCORES

log "Validating sample metadata"
# Fails early on malformed metadata (duplicate IDs, bad headers) before hours of compute.
have "$D/metadata.qzv" || run qiime metadata tabulate \
    --m-input-file "$METADATA" \
    --o-visualization "$D/metadata.qzv"

if [[ "$RUN_FASTQC" == "yes" ]]; then
    if command -v fastqc >/dev/null && command -v multiqc >/dev/null; then
        log "FastQC + MultiQC on raw reads (outlier samples, per-base quality)"
        mkdir -p "$D/fastqc"
        have "$D/multiqc_report.html" || {
            run fastqc -t "$NCORES" -o "$D/fastqc" "$RAW_DIR"/*.f*q.gz
            run multiqc "$D/fastqc" --outdir "$D" --filename multiqc_report.html
        }
    else
        log "fastqc/multiqc not on PATH; skipping (install them in a separate env)."
    fi
fi

log "Importing reads (${PLATFORM}, ${LAYOUT}, ${IMPORT_MODE})"
if [[ "$PLATFORM" == "pacbio" || "$LAYOUT" == "single" ]]; then
    TYPE='SampleData[SequencesWithQuality]'
    MANIFEST_FMT=SingleEndFastqManifestPhred33V2
else
    TYPE='SampleData[PairedEndSequencesWithQuality]'
    MANIFEST_FMT=PairedEndFastqManifestPhred33V2
fi

if ! have "$D/reads.qza"; then
    case "$IMPORT_MODE" in
        casava)
            run qiime tools import \
                --type "$TYPE" \
                --input-path "$RAW_DIR" \
                --input-format CasavaOneEightSingleLanePerSampleDirFmt \
                --output-path "$D/reads.qza" ;;
        manifest)
            [[ -f "$MANIFEST" ]] || die "Manifest not found: $MANIFEST"
            run qiime tools import \
                --type "$TYPE" \
                --input-path "$MANIFEST" \
                --input-format "$MANIFEST_FMT" \
                --output-path "$D/reads.qza" ;;
        *) die "IMPORT_MODE must be casava or manifest" ;;
    esac
fi

have "$D/reads_summary.qzv" || run qiime demux summarize \
    --i-data "$D/reads.qza" \
    --o-visualization "$D/reads_summary.qzv"

log "CHECKPOINT: open $D/reads_summary.qzv (and multiqc_report.html)."
log "Confirm sample count and per-sample depth; note where quality drops."
