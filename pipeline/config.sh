# shellcheck shell=bash disable=SC2034  # variables are consumed by the sourcing stage scripts
# -----------------------------------------------------------------------------
# Project configuration for the combined amplicon SOP.
#
# Every stage script sources this file. Values left empty ("") are decisions
# that must be made after looking at a visualization produced by an earlier
# stage; the stage that needs the value stops and names the .qzv to inspect.
#
# Source labels in brackets refer to docs/01-workflow-comparison.md, section 1.
#
# Every stage prints the settings it used at the start of its log, each tagged
# (config), (preset ...) or (yours; preset ... would set ...), and saves the same
# list to results/<stage>/settings_used.tsv.
# -----------------------------------------------------------------------------

# ---- Preset -----------------------------------------------------------------
# ""        : the generic defaults below (written for Illumina 16S)
# "18S_PR2" : 18S V4 + PR2 defaults; see the preset block at the end of this file.
#             A preset only changes settings you left at the value shipped here;
#             anything you edited yourself is kept. To keep a value that equals the
#             shipped default, use the "Your overrides" section at the end.
PRESET=""

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
# Rank names of your reference taxonomy, in order. This is the single place ranks are
# defined: QIIME 2 collapses by number, R looks the name up here (via the stage-10 export).
TAXONOMY_RANKS="Domain Phylum Class Order Family Genus Species"   # SILVA / GTDB / Greengenes2
OVERVIEW_RANK_LEVEL=2                     # rank for overview barplots (2 = Phylum here)
# Stage 05a can train a PR2 classifier with RESCRIPt, which also guarantees the
# scikit-learn version matches this QIIME 2 release.
BUILD_PR2_CLASSIFIER="no"                 # yes = run 05a_build_pr2_classifier.sh
PR2_VERSION="5.1.1"                       # 5.1.1 | 5.1.0 | 5.0.0 (rescript get-pr2-data)
# A pre-trained classifier MUST match the scikit-learn version of your QIIME 2
# release (see SOP, stage 5).
CLASSIFIER=""                             # path to a classifier .qza
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
# Rank number (1-based, in TAXONOMY_RANKS) for collapsed-table differential abundance,
# used by QIIME 2 (taxa collapse) and R alike. 6 = Genus for the 7 ranks above.
TAXA_COLLAPSE_LEVEL=6

# =============================================================================
# Preset block. Format:  preset VARIABLE "value shipped above" "preset value"
# The middle argument must match the value shipped above (tests/smoke_test.sh checks
# this); if you changed the variable, your value is kept and the logs say so.
# =============================================================================
if [[ "$PRESET" == "18S_PR2" ]]; then
    preset MARKER               "16S"   "18S"
    preset PRIMER_F             "GTGYCAGCMGCCGCGGTAA"  "CYGCGGTAATTCCAGCTC"     # 18S V4
    preset PRIMER_R             "CCGYCAATTYMTTTRAGTTT" "AYGGTATCTRATCRTCTTYG"
    preset DADA2_POOLING        "independent" "pseudo"   # helps rare ASVs seen in several samples
    preset TREE_METHOD          "auto"  "none"            # no reliable tree across all eukaryotes
    preset RUN_KMER_DIVERSITY   "no"    "yes"             # tree-free relatedness-aware metrics
    preset MIN_SAMPLES          "2"     "1"               # keep taxa found at a single site
    # Rare-ASV frequency rule: keep ONE of the next two lines active.
    preset MIN_FREQ_FRACTION    "0.001" "0.001"           # keep the 0.1% rule
    # preset MIN_FREQ_FRACTION  "0.001" "0"               # drop the 0.1% rule
    preset TAXONOMY_RANKS       "Domain Phylum Class Order Family Genus Species" \
                                "Domain Supergroup Division Subdivision Class Order Family Genus Species"
    preset OVERVIEW_RANK_LEVEL  "2"     "3"               # Division
    preset TAXA_COLLAPSE_LEVEL  "6"     "8"               # Genus (9 = Species, 7 = Family, ...)
    preset BUILD_PR2_CLASSIFIER "no"    "yes"
    preset CLASSIFIER           ""      "${OUT}/05a_pr2_classifier/pr2_${PR2_VERSION}_classifier.qza"
    preset CLASSIFIER_LABEL     "reference" "PR2_${PR2_VERSION}"
    # Prokaryotes are matched on the domain label so the diatom genus Bacteriastrum is kept;
    # organelle tags are PR2's (:plas plastid, :mito, :apic apicoplast, :chro chromatophore).
    # PR2 also tags nucleomorphs (:nucl); add ",:nucl" to drop them too.
    preset TAXA_EXCLUDE         "mitochondria,chloroplast" \
                                "d__Bacteria,d__Archaea,:plas,:mito,:apic,:chro"
fi

# =============================================================================
# Your overrides. Format:  override VARIABLE "value"
# Applied last, so these always win over the preset and the values above, even when
# your value equals a default (e.g. keeping MIN_SAMPLES=2 with PRESET="18S_PR2").
# Each stage reports them as (yours: override ...). Example:
#   override MIN_SAMPLES "2"
# =============================================================================
