# shellcheck shell=bash disable=SC2034  # variables are consumed by the sourcing stage scripts
# -----------------------------------------------------------------------------
# Project configuration for the combined amplicon SOP.
#
# Every stage script sources this file. Values left empty ("") are decisions
# that must be made after looking at a visualization produced by an earlier
# stage; the stage that needs the value stops and names the .qzv to inspect.
#
# Source labels in brackets refer to docs/01-workflow-comparison.md, section 1.
# -----------------------------------------------------------------------------

# ---- Environment ------------------------------------------------------------
# QIIME 2 release this project runs on. Supported here: 2026.1, 2026.4, 2026.7.
# It selects version-specific behaviour (ANCOM-BC2 visualizer, q2-boots seed).
Q2_EPOCH="2026.7"
NCORES=4                                  # threads for parallel steps [MH_*]
RANDOM_SEED=711                           # resampling seed (2026.4+) [MH_2026.3-R uses 711]

# ---- Paths ------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${PROJECT_DIR}/raw_data"         # demultiplexed FASTQ(.gz) files
METADATA="${PROJECT_DIR}/metadata.tsv"    # QIIME 2 metadata (first column = sample id)
OUT="${PROJECT_DIR}/results"              # all outputs land here, one folder per stage

# ---- Study design -----------------------------------------------------------
PLATFORM="illumina"                       # illumina | pacbio
LAYOUT="paired"                           # paired | single (Illumina only; PacBio is always single)
MARKER="16S"                              # 16S | 18S | ITS | other
GROUP_COL=""                              # main categorical metadata column, e.g. "treatment"
ADONIS_FORMULA=""                         # e.g. "treatment+site" (terms are fitted sequentially)
DA_REFERENCE=""                           # e.g. "treatment::control" for ANCOM-BC2 reference level
DA_SAMPLE_FILTER=""                       # optional SQL where-clause, e.g. "[treatment] IN ('A', 'B')"

# ---- 01 Import --------------------------------------------------------------
# casava   : files named like 105CHE6WT_S325_L001_R1_001.fastq.gz      [MH_qiime2-amplicon-2024.5]
# manifest : any naming; supply a QIIME 2 fastq manifest (TSV)          [experienced_researchers]
IMPORT_MODE="casava"
MANIFEST="${PROJECT_DIR}/manifest.tsv"
RUN_FASTQC="yes"                          # needs fastqc + multiqc on PATH (not part of QIIME 2)

# ---- 02 Primer removal (Illumina) ------------------------------------------
# Leave both empty to skip trimming (e.g. EMP protocol reads start after the primer).
# Common primer sets from [MH_2026.3-workflow]:
#   16S V4-V5 (universal) : GTGYCAGCMGCCGCGGTAA / CCGYCAATTYMTTTRAGTTT
#   16S V6-V8 (bacteria)  : ACGCGHNRAACCTTACC   / ACGGGCRGTGWGTRCAA
#   16S V6-V8 (archaea)   : TYAATYGGANTCAACRCC  / CRGTGWGTRCAAGGRGCA
#   18S V4                : CYGCGGTAATTCCAGCTC  / AYGGTATCTRATCRTCTTYG
#   ITS2 (fungi)          : GTGAATCATCGAATCTTTGAA / TCCTCCGCTTATTGATATGC
#   16S V4 (EMP 515F/806R): GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT
PRIMER_F="GTGYCAGCMGCCGCGGTAA"
PRIMER_R="CCGYCAATTYMTTTRAGTTT"

# ---- 03 Denoising -----------------------------------------------------------
DENOISER="dada2"                          # dada2 (default) | deblur (Illumina 16S/18S/ITS)
# DADA2 (Illumina). Choose from the quality plot of the trimmed reads; make sure
# trunc_f + trunc_r >= amplicon length (after primer removal) + 12 nt overlap.
DADA2_TRIM_LEFT_F=0
DADA2_TRIM_LEFT_R=0
DADA2_TRUNC_F=""                          # e.g. 240  (single-end: only this one is used)
DADA2_TRUNC_R=""                          # e.g. 200  (use 0 for both on ITS: variable length)
DADA2_POOLING="independent"               # independent | pseudo
# Deblur: trim length from the length-percentile table of the merged reads.
DEBLUR_TRIM_LENGTH=""                     # e.g. 366 for V4-V5 [MH_2026.3-workflow]
DEBLUR_REF_SEQS=""                        # required for 18S/ITS: positive-filter reference (.qza)
# PacBio HiFi full-length (dada2 denoise-ccs) [MH_2026.3-workflow]
CCS_FRONT="AGRGTTYGATYMTGGCTCAG"          # 27F (bacteria); see SOP for archaea/18S/ITS sets
CCS_ADAPTER="RGYTACCTTGTTACGACTT"         # 1492R
CCS_MIN_LEN=1200
CCS_MAX_LEN=1800
CCS_MAX_EE=18

# ---- 04 Feature filtering ---------------------------------------------------
MIN_SAMPLES=2                             # keep ASVs seen in >= N samples [gut_to_soil_tutorial]
MIN_FREQ_FRACTION="0.001"                 # keep ASVs with total count >= fraction x mean depth [MH_*]
                                          # set to 0 to disable the frequency rule

# ---- 05 Taxonomy ------------------------------------------------------------
# MUST match the scikit-learn version of your QIIME 2 release (see SOP, stage 5).
CLASSIFIER=""                             # path to a pre-trained classifier .qza
CLASSIFIER_LABEL="reference"              # short name shown in tabulate-seqs, e.g. "SILVA138"
CLASSIFY_CONFIDENCE="0.7"

# ---- 06 Taxonomy / sample filtering -----------------------------------------
TAXA_EXCLUDE="mitochondria,chloroplast"   # 16S default [MH_*]; "" to skip
TAXA_INCLUDE=""                           # e.g. "p__" to drop phylum-unassigned ASVs (not default)
MIN_SAMPLE_DEPTH=""                       # set after viewing rarefaction curves (e.g. 1000-4000)

# ---- 07 Phylogeny -----------------------------------------------------------
# auto: 16S -> sepp ; 18S -> denovo ; ITS/other -> none
TREE_METHOD="auto"                        # auto | sepp | denovo | none
SEPP_REF=""                               # SEPP reference .qza (e.g. GG2 2022.10 backbone) [MH_2026.3-workflow]

# ---- 08 Diversity -----------------------------------------------------------
SAMPLING_DEPTH=""                         # even-sampling depth; set after stage 06 review
BOOTS_N=100                               # resampling iterations [gut_to_soil_tutorial (dev)]
RUN_KMER_DIVERSITY="no"                   # yes = also run tree-free kmer-diversity [gut_to_soil_tutorial]

# ---- 09 Composition / differential abundance --------------------------------
TAXA_COLLAPSE_LEVEL=6                     # 6 = genus for 7-rank taxonomies
