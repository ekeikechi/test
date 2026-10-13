# Combined amplicon SOP (QIIME 2 2026.7 + R)

This SOP merges the eight sources compared in
[`01-workflow-comparison.md`](01-workflow-comparison.md) into one annotated workflow. It
starts from raw demultiplexed reads and ends with figures and statistics. Each step cites
the source it came from using the labels you gave (`[MH_2026.3-workflow]`, `[gut_to_soil_tutorial]`, …).
Where the SOP departs from every source, it says so and gives the reason.

The commands are implemented in [`../pipeline/`](../pipeline) (ten stage scripts plus
`config.sh`) and [`../R/downstream_analysis.R`](../R/downstream_analysis.R). This document
explains the decisions; the scripts carry them out.

---

## 0. Design decisions at a glance

| # | Decision | Taken from | Why |
|---|---|---|---|
| D1 | Target **QIIME 2 2026.7**; notes for 2026.1 | `gut_to_soil_tutorial` (2026.7), `MH_2026.3-workflow` (2026.1) | Current stable release. 2026.1 differences are handled by `Q2_EPOCH`. |
| D2 | Import with Casava naming **or** a manifest | `MH_2026.3-workflow`, `experienced_researchers` | The manifest "will always work". |
| D3 | Raw QC: FastQC + MultiQC **and** `demux summarize` | `MH_*`, all QIIME 2 tutorials | MultiQC catches outlier samples; `demux summarize` drives the truncation choice. |
| D4 | Primer removal with `cutadapt --p-discard-untrimmed --p-no-indels` (Illumina) | `MH_*` | Removes off-target reads and primer dimers. Skipped when primers were not sequenced (EMP protocol). |
| D5 | **DADA2** is the default denoiser; Deblur is an option; DADA2 `denoise-ccs` for PacBio | `gut_to_soil_tutorial`, `experienced_researchers`; `MH_*` for Deblur and CCS | DADA2 needs no merge or q-filter step, avoids losing reads to stringent merging, and works for every marker. |
| D6 | Rare-ASV filter (≥2 samples **and** ≥0.1% of mean depth) **before** taxonomy | `gut_to_soil_tutorial` + `MH_*`, `conceptual_overview` | One documented rule; faster classification. |
| D7 | Version-matched, full-length NB classifier (environment-weighted when one fits); `--p-read-orientation same` | `MH_2026.3-workflow`, `conceptual_overview`, QIIME 2 resources page | Avoids sklearn version errors; weighted classifiers improve accuracy. |
| D8 | Remove mitochondria/chloroplast; keep unclassified ASVs | `MH_*` | Unclassified ASVs are diagnostic, not noise by default. |
| D9 | Remove failed samples using rarefaction curves; check losses per group | `MH_*`, `gut_to_soil_tutorial` | Explicit, reportable threshold. |
| D10 | 16S: SEPP **+ `fragment-insertion filter-features`**; 18S: MAFFT/FastTree pipeline; ITS: no tree | `MH_2026.3-workflow`, `moving_pictures_tutorial`, `gut_to_soil_tutorial` | The filter step is missing in MH (comparison E7). |
| D11 | Diversity via **q2-boots core-metrics** (rarefaction, n = 100, fixed seed) | `gut_to_soil_tutorial`, `conceptual_overview` | Averages over many rarefactions instead of one random draw. |
| D12 | PERMANOVA (pairwise) **+ PERMDISP** + multi-factor `adonis` | `moving_pictures_tutorial`, `MH_2026.3-stats-visualization`, `MH_2026.3-R` | Not in any source: the dispersion check is needed to interpret PERMANOVA. |
| D13 | ANCOM-BC2 in QIIME 2 (ASV and genus); R consensus of ANCOM-BC2 + ALDEx2 + MaAsLin3 | `gut_to_soil_tutorial`, `MH_2026.3-R` | Current QIIME 2 method plus MH's multi-tool practice, with the same input to every tool. |
| D14 | Export for R plus `replay-provenance` and `replay-supplement` | `MH_*`, `gut_to_soil_tutorial` | Portable data plus an exact methods record. |

---

## 1. Installation

**QIIME 2 2026.7** (Linux). On macOS, add `CONDA_SUBDIR=osx-64` in front of the create
command; there is no osx-arm64 build for this release.

```bash
conda env create --name rachis-qiime2-2026.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml
conda activate rachis-qiime2-2026.7
qiime info
```

For **2026.1** (the release MH 2026.3 was written for), use
`…/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml` and name the environment
`qiime2-amplicon-2026.1`. Then set `Q2_EPOCH="2026.1"` in `pipeline/config.sh`.

**FastQC/MultiQC** are not part of QIIME 2 `[MH_*]`. Install them in their own environment
and make sure they are on `PATH` when stage 01 runs:

```bash
conda create -n readqc -c conda-forge -c bioconda fastqc multiqc
```

**R** (stage 11): R ≥ 4.3 with `phyloseq`, `vegan`, `ggplot2`, `dplyr`, `tidyr`, `ALDEx2`,
`maaslin3`, `ANCOMBC`, `ape`. For example:

```bash
conda create -n amplicon-r -c conda-forge -c bioconda r-base=4.5 bioconductor-phyloseq r-vegan \
  r-ggplot2 r-dplyr r-tidyr bioconductor-aldex2 bioconductor-maaslin3 bioconductor-ancombc r-ape
```

## 2. Project layout and how to run

```
project/
├── raw_data/            demultiplexed FASTQ(.gz)
├── metadata.tsv         QIIME 2 metadata (first column = sample id; no personal data)
├── manifest.tsv         only when IMPORT_MODE=manifest
├── pipeline/            config.sh, lib.sh, stage scripts 01_…10_ (+ optional 05a_)
├── R/downstream_analysis.R
└── results/             created by the pipeline: one folder per stage
```

Edit `pipeline/config.sh`, then run the stages in order:

```bash
bash pipeline/01_import_qc.sh
bash pipeline/02_trim_primers.sh
# ...03, 04, 05a (only builds a PR2 classifier when BUILD_PR2_CLASSIFIER=yes), 05...
bash pipeline/10_export_provenance.sh
Rscript R/downstream_analysis.R results/10_export results/11_R   # group from GROUP_COL
```

The pipeline needs human decisions at several points: truncation lengths, the minimum sample
depth, and the even-sampling depth. When a decision variable is still empty in `config.sh`,
the stage that needs it **stops and names the `.qzv` to inspect**. Set the value and re-run
that stage. Stages skip outputs that already exist, so re-running is safe. To redo a step,
delete its output files yourself first.

View `.qzv` files with QIIME 2 View (<https://view.qiime2.org>) `[all]`.

### 2.1 Presets, the settings report, and one place for shared settings

- **`PRESET="18S_PR2"`** (top of `config.sh`) switches the generic Illumina-16S defaults
  to 18S V4 + PR2. The preset block at the end of `config.sh` lists every value it changes.
  A preset **only changes a setting you left at the value shipped in `config.sh`**: if you
  edited it, your value is kept. A preset cannot tell a deliberate choice from an untouched
  default, so to keep a value that *equals* the shipped default (e.g. `MIN_SAMPLES=2` with
  18S_PR2), write `override MIN_SAMPLES "2"` in the **Your overrides** section at the end of
  `config.sh`; overrides are applied last and always win. The rare-ASV frequency rule has two lines in the block;
  keep one active (0.001 keeps the 0.1% rule, which is active as shipped; 0 drops it).
- **Settings report.** Each stage starts its log with the settings it used, tagged
  `(config)`, `(preset 18S_PR2)`, `(yours; preset 18S_PR2 would set '…')` or
  `(yours: override …)`. It saves the
  same list to `results/<stage>/settings_used.tsv`. The R script prints and saves its own
  list (`settings_used.csv`).
- **Shared settings live only in `config.sh`.** These are `TAXONOMY_RANKS` (rank names),
  `TAXA_COLLAPSE_LEVEL` (rank number for collapsed-table differential abundance),
  `OVERVIEW_RANK_LEVEL`, `SAMPLING_DEPTH`, `GROUP_COL`, `DA_REFERENCE` and `RANDOM_SEED`.
  Stage 10 writes them to `results/10_export/analysis_settings.tsv`, which the R script
  reads, so QIIME 2 and R always analyse the same rank and depth.
- **Safety checks.** A collapse level beyond the number of ranks stops every stage.
  Stage 10 and the R script stop if the taxonomy has more ranks than `TAXONOMY_RANKS`
  instead of silently shifting rank labels.
- **Old bash.** The scripts run on bash 3.2+ (tested on 4.3.48 and 5.2). Arrays that may be
  empty are expanded so that bash < 4.4 does not fail with "unbound variable".

---

## Stage 1. Metadata, raw-read QC and import (`01_import_qc.sh`)

1. **Validate metadata** with `qiime metadata tabulate` before anything else
   `[gut_to_soil_tutorial]` `[moving_pictures_tutorial]`. It catches duplicate IDs and bad
   headers before hours of computing. Never put personal identifiers in metadata: QIIME 2
   stores metadata in provenance `[moving_pictures_tutorial]`.
2. **FastQC + MultiQC** on the raw FASTQs `[MH_qiime2-amplicon-2024.5 §1.3]`. Look for outlier
   samples and check that per-base quality is mostly above Q30.
3. **Import**:
   - `IMPORT_MODE=casava` for Illumina names like `105CHE6WT_S325_L001_R1_001.fastq.gz`
     (`CasavaOneEightSingleLanePerSampleDirFmt`) `[MH_*]`.
   - `IMPORT_MODE=manifest` for anything else, including PacBio and SRA downloads
     (`PairedEndFastqManifestPhred33V2` / `SingleEndFastqManifestPhred33V2`)
     `[MH_2026.3-workflow]` `[experienced_researchers]`. A manifest is a TSV with columns
     `sample-id`, `forward-absolute-filepath`, and (for paired data) `reverse-absolute-filepath`.
   - Multiplexed data must be demultiplexed first (`demux emp-*` or `cutadapt demux-*`).
     Dual-index paired-end data must be demultiplexed outside QIIME 2
     `[experienced_researchers]`. Ask the sequencing centre for demultiplexed data
     `[gut_to_soil_tutorial]`.
4. **`demux summarize`** on the imported reads `[MH_2026.3-workflow §2.2]`.

**Checkpoint:** sample count and per-sample depth look as expected.

## Stage 2. Primer removal (`02_trim_primers.sh`, Illumina)

`cutadapt trim-paired` (or `trim-single`) with `--p-discard-untrimmed --p-no-indels`
`[MH_*]`. Primer sets in `config.sh` `[MH_2026.3-workflow]`:

| Target | Forward | Reverse |
|---|---|---|
| 16S V4–V5 ("universal") | `GTGYCAGCMGCCGCGGTAA` | `CCGYCAATTYMTTTRAGTTT` |
| 16S V6–V8 (bacteria) | `ACGCGHNRAACCTTACC` | `ACGGGCRGTGWGTRCAA` |
| 16S V6–V8 (archaea) | `TYAATYGGANTCAACRCC` | `CRGTGWGTRCAAGGRGCA` |
| 18S V4 | `CYGCGGTAATTCCAGCTC` | `AYGGTATCTRATCRTCTTYG` |
| ITS2 (fungi) | `GTGAATCATCGAATCTTTGAA` | `TCCTCCGCTTATTGATATGC` |

If the protocol does not sequence the primers (EMP 515F/806R as used in both QIIME 2
tutorials), leave `PRIMER_F`/`PRIMER_R` empty and the stage passes the reads through
unchanged. PacBio primers are removed inside `denoise-ccs` (stage 3) `[MH_2026.3-workflow]`.

**Checkpoint:** most reads survive trimming. A large loss almost always means a wrong
primer sequence `[MH_2026.3-workflow §2.3.2]`.

## Stage 3. Denoising (`03_denoise.sh`)

### 3A. Illumina, DADA2 (default)

`dada2 denoise-paired` merges pairs, removes chimeras and filters phiX internally. **Do
not** merge or q-score-filter first `[experienced_researchers]` `[gut_to_soil_tutorial]`.

Choose truncation positions in the *Interactive Quality Plot* of
`02_trim/reads_trimmed_summary.qzv`. Truncate where median quality drops
`[moving_pictures_tutorial]` `[gut_to_soil_tutorial]`. Reads are dropped if they are
shorter than the truncation length, and pairs need **≥12 nt overlap** (the DADA2 default
`min_overlap = 12`), so:

```
trunc_len_f + trunc_len_r  ≥  amplicon length after primer removal + 12
```

Example for 16S V4–V5 (515F/926R): the amplicon is about 411 bp, or about 372 nt after
primer removal (MH reports a median merged length of 372 nt). So `240 + 160 = 400 ≥ 384`
works.

For **18S V4**, amplicon length varies across eukaryote groups. If the truncation lengths
are too short for the longer amplicons, those reads fail to merge and their taxa vanish
without an error. Check the *merged* column per sample in `denoising_stats.qzv`.

For **ITS**, amplicon length varies a lot. Use `trunc_len 0` for both reads and rely on
`max_ee` (general DADA2 practice for ITS, not taken from the sources).

### 3B. Illumina, Deblur (option for 16S, or when speed matters)

`vsearch merge-pairs` → `quality-filter q-score` → `deblur denoise-16S`
`[MH_qiime2-amplicon-2024.5 §2]` `[MH_2026.3-workflow §2.3.3]` `[moving_pictures_tutorial]`.

- Pick `--p-trim-length` from the length-percentile table in `reads_qfilt_summary.qzv`.
  Take the shortest length that keeps ~98% of reads, if it is plausible for the amplicon
  (e.g. 366 nt for V4–V5) `[MH_2026.3-workflow]`.
- `--p-min-reads 1` turns off Deblur's across-sample abundance filter (default 10), so that
  stage 4 applies the same rare-ASV rule whichever denoiser was used. Per-sample singletons
  are still removed by `--p-min-size` (default 2) (comparison E6).
- 18S/ITS need `denoise-other` with a positive-filter reference (PR2 for 18S, UNITE for ITS)
  `[MH_*]`. `denoise-16S` positive-filters against Greengenes 13_8 85% OTUs
  `[experienced_researchers]`, which is one reason DADA2 is the default here.
- If few reads merge, `vsearch merge-pairs` is stricter than PEAR `[MH_*]`. Switch to DADA2
  rather than exporting and re-importing.

### 3C. PacBio HiFi full-length, DADA2 `denoise-ccs` `[MH_2026.3-workflow]`

Primer removal, re-orientation and length gating all happen inside DADA2. No quality
filter is applied (HiFi reads are ≥99% accurate). Parameter sets:

| Target | `--p-front` | `--p-adapter` | min–max length | `--p-max-ee` |
|---|---|---|---|---|
| 16S bacteria (27F/1492R) | `AGRGTTYGATYMTGGCTCAG` | `RGYTACCTTGTTACGACTT` | 1200–1800 | 18 |
| 16S archaea (in testing at MH) | `TCCGGTTGATCCYGCCGG` | `CRGTGWGTRCAAGGRGCA` | 1200–1800 | 18 |
| 18S full-length | `CTGGTTGATYCTGCCAGT` | `TGATCCTTCTGCAGGTTCACCTAC` | 1000–3000 | 30 |
| ITS full-length (fungi) | `TAGAGGAAGTAAAAGTCGTAA` | `TCCTCCGCTTWTTGWTWTGC` | 300–900 | 9 |

All use `--p-pooling-method pseudo --p-n-reads-learn 100000`. The 16S window starts at
1200 nt as in MH's command; MH's text says 1300 (comparison E5).

**Checkpoint (all branches):** in `denoising_stats.qzv`, which the pipeline merges with the
sample metadata `[gut_to_soil_tutorial]`, check the fraction of reads retained.
- Losses at *filtered*: truncation is too aggressive.
- Losses at *merged*: overlap is too short.
- Losses at *non-chimeric*: often leftover primers.

For Deblur, check `deblur_stats.qzv`. Very low post-Deblur depth is a red flag `[MH_*]`.

## Stage 4. Rare-ASV filtering (`04_filter_features.sh`)

Before classification, which makes classification faster and lighter on memory
`[gut_to_soil_tutorial]` `[conceptual_overview]`:

- **Prevalence**: keep ASVs found in ≥ `MIN_SAMPLES` (default 2) samples
  `[gut_to_soil_tutorial]`.
- **Frequency**: keep ASVs with a total count ≥ round(0.001 × mean sample depth). This is the
  0.1% Illumina bleed-through rule `[MH_*]`. The script computes the mean from the table
  and records the threshold in `filter_thresholds.tsv` for the methods section.

Set `MIN_SAMPLES=1` for designs with unique, unreplicated samples. Set
`MIN_FREQ_FRACTION=0` to drop the frequency rule. Don't add a second, different abundance
filter in R (MH-R's ≥100 reads): one documented rule is easier to defend.

**Checkpoint:** many ASVs but few reads should be removed.

## Stage 5. Taxonomy (`05_taxonomy.sh`)

`classify-sklearn` with a **pre-trained, full-length** Naive Bayes classifier
`[MH_2026.3-workflow]` `[conceptual_overview]`. Region-specific classifiers are no longer
recommended by MH or by QIIME 2.

**The classifier must match your QIIME 2 release's scikit-learn version.** Otherwise
`classify-sklearn` stops with a `ValueError` (verified in q2-feature-classifier source):

| QIIME 2 | scikit-learn | Where to get classifiers |
|---|---|---|
| 2024.5 – 2026.1 | 1.4.2 | QIIME 2 builds linked from `MH_2026.3-workflow`: `https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza`, `.../greengenes2/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza`, `.../gtdb/gtdb_classifier_r220.qza`, plus diverse- and human-stool-**weighted** SILVA/GTDB builds |
| 2026.4 – 2026.7 | 1.7.1 | QIIME 2 now points to external hosts: **SILVA** at <https://www.arb-silva.de/documentation/classifiers/qiime-2>, **GTDB** at Zenodo [10.5281/zenodo.19446101](https://doi.org/10.5281/zenodo.19446101). Or train your own with RESCRIPt `[gut_to_soil_tutorial]`. |

- **PR2 (18S), stage 05a.** With `BUILD_PR2_CLASSIFIER=yes` (set by the 18S_PR2 preset),
  `05a_build_pr2_classifier.sh` trains a classifier in your own QIIME 2 environment, so the
  scikit-learn version always matches. It runs `rescript get-pr2-data --p-version 5.1.1`
  → `rescript dereplicate --p-mode uniq --p-rank-handles disable` →
  `feature-classifier fit-classifier-naive-bayes`, on the full-length reference.
  `--p-rank-handles disable` matters: RESCRIPt otherwise backfills missing ranks with its
  default 7 NCBI-style prefixes. Training on all of PR2 needs a large-memory node.
  - RESCRIPt labels the 9 PR2 ranks `d__ sgr__ dv__ dvs__ c__ o__ f__ g__ s__`; set
    `TAXONOMY_RANKS` to the 9 names (the preset does). PR2 fills unknown ranks with
    placeholders such as `g__Trichostomatia_XX`, so some "genera" are not real genera.
  - Choose `TAXA_COLLAPSE_LEVEL` by rank number: 1 domain, 2 supergroup, 3 division,
    4 subdivision, 5 class, 6 order, 7 family, **8 genus**, 9 species. Short V4 reads often
    fail the 0.7 confidence at species level, so check in `rep_seqs_taxonomy.qzv` how many
    reads reach the rank you choose.
- Prefer an **environment-weighted** classifier when one matches your sample type
  `[conceptual_overview]` `[moving_pictures_tutorial]`.
- Other references: UNITE (ITS), PR2 (18S), HOMD (oral), MitoFish (12S) `[MH_2026.3-workflow]`.
  `MH_2026.3-workflow` shows how to build a HOMD classifier with `fit-classifier-naive-bayes`.
- Classifiers contain executable code. Use only ones from sources you trust
  `[gut_to_soil_tutorial (dev)]`.
- Memory: 8–16 GB or more `[conceptual_overview]`. If you run out, lower
  `--p-reads-per-batch` and `--p-n-jobs` `[MH_*]`.
- `--p-read-orientation same` is set because the reads are already in amplicon orientation
  after primer trimming or `denoise-ccs` `[MH_2026.3-workflow]`.

**Checkpoint:** open `rep_seqs_taxonomy.qzv`, which shows sequence, taxonomy, frequency
and BLAST links `[gut_to_soil_tutorial]`. BLAST about 5 abundant ASVs from different phyla
and compare `[MH_*]` `[moving_pictures_tutorial]`. Trust the curated classifier over
individual BLAST hits `[gut_to_soil_tutorial]`. Many phylum-unassigned ASVs mean a wrong
classifier for the marker, wrong orientation, or off-target amplification.

## Stage 6. Taxonomy filtering and failed samples (`06_filter_taxa_samples.sh`)

1. `taxa filter-table --p-exclude mitochondria,chloroplast` for 16S `[MH_*]`. Matching is a
   case-insensitive substring (verified). **Keep unclassified ASVs by default**: a high
   fraction is informative `[MH_*]`. Check the labels your reference uses. SILVA uses
   *Mitochondria* (family) and *Chloroplast* (order). GTDB has no organelle sequences, so
   the filter matches nothing and organelle reads may be classified elsewhere.
   **PR2 / 18S** (preset): `d__Bacteria,d__Archaea,:plas,:mito,:apic,:chro`. PR2 5.1.1 tags
   organelle sequences with `:plas`, `:mito`, `:apic`, `:chro`, plus `:nucl` (nucleomorph),
   which you can add; this was checked in the 5.1.1 taxonomy file. Match prokaryotes on the
   domain label: a bare `Bacteria` would also remove the diatom genus *Bacteriastrum*. For
   ITS, set `TAXA_EXCLUDE=""` or exclude host lineages.
2. **Rarefaction curves** up to the median depth `[moving_pictures_tutorial]`
   `[gut_to_soil_tutorial]`, drawn twice: grouped by metadata and per sample
   `[MH_qiime2-amplicon-2024.5 §6]`. MH notes that QIIME 2 shows only grouped curves when
   metadata is given, hence the second, per-sample run.
3. Set `MIN_SAMPLE_DEPTH` to remove **failed** samples. Researchers often use 1000–4000, or
   much lower to drop only total failures `[MH_*]`. Check how many samples each group loses
   (*Interactive Sample Detail*, grouped by your main column) `[gut_to_soil_tutorial]`.
   Use `1` to keep everything.
4. Filter the representative sequences to the final table.

**Checkpoint:** choose `SAMPLING_DEPTH` (the even-sampling depth) in
`table_final_summary.qzv`. Make it as high as possible while losing few samples, near where
the curves level off `[all]`. Results are fairly stable across reasonable depths, so try
two and compare `[gut_to_soil_tutorial (dev)]`.

## Stage 7. Phylogeny (`07_phylogeny.sh`)

| Marker | Method | Source |
|---|---|---|
| 16S | `fragment-insertion sepp` into a reference tree, **then `fragment-insertion filter-features`** to drop ASVs SEPP could not place | `[MH_2026.3-workflow §5.1]`, `[experienced_researchers]` (reference trees are higher quality) |
| 18S | `phylogeny align-to-tree-mafft-fasttree` (one command doing MAFFT → mask → FastTree → midpoint root) | `[moving_pictures_tutorial]`; same steps as `[MH_*]` §5B |
| ITS / other | no tree: non-phylogenetic metrics, optionally kmer-based diversity | `[gut_to_soil_tutorial]` |

SEPP references `[MH_2026.3-workflow]`: Greengenes2 2022.10 backbone (MH's current choice,
from the MH server), `https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-silva-128.qza`,
or `.../sepp-refs-gg-13-8.qza`.

## Stage 8. Diversity and group statistics (`08_diversity.sh`)

**`boots core-metrics`** rarefies `BOOTS_N` (default 100) times without replacement. For
each metric it then takes the median alpha value and the medoid distance matrix, and builds
Emperor plots and a combined scatter plot `[gut_to_soil_tutorial]` `[conceptual_overview]`.
This replaces the single-rarefaction `core-metrics-phylogenetic` used by `[MH_*]` and
`[moving_pictures_tutorial]`.

Metrics: Shannon, observed features, Pielou evenness, Bray–Curtis and Jaccard. With a tree,
it adds Faith PD plus weighted and unweighted UniFrac. On 2026.4+ the run is reproducible
via `--p-random-seed`.

The stable gut-to-soil page used bootstrapping (`replacement=True`); its development
version and this SOP use rarefaction, which is the more widely used approach. Optionally
(`RUN_KMER_DIVERSITY=yes`), `boots kmer-diversity` gives tree-free, relatedness-aware
metrics `[gut_to_soil_tutorial]`. This is useful for ITS. Remember that "observed features"
then counts kmers, not ASVs.

Statistics:

- `diversity alpha-group-significance` (Kruskal–Wallis) for categorical columns and
  `alpha-correlation` (Spearman) for numeric ones `[moving_pictures_tutorial]`
  `[MH_2026.3-stats-visualization]`. The q2-stats equivalent is deprecated in 2026.7.
- `diversity beta-group-significance --p-pairwise` with **PERMANOVA and PERMDISP** for
  `GROUP_COL` `[moving_pictures_tutorial]` `[MH_2026.3-stats-visualization]`. PERMDISP is an
  SOP addition: a significant PERMANOVA together with a significant PERMDISP may reflect a
  difference in spread rather than in centroid.
- `diversity adonis` for multi-factor models (`ADONIS_FORMULA`). QIIME 2 fits terms
  **sequentially** (`by = "terms"`, so order matters); the R script reports **marginal**
  terms (comparison E14).

## Stage 9. Composition and differential abundance (`09_composition_da.sh`)

- `taxa barplot` on the ASV table `[all]`. Also a per-group summed barplot through
  `feature-table group --p-mode sum` `[MH_qiime2-amplicon-2024.5 §7]`; use it for
  overview figures only.
- **ANCOM-BC2** on **unrarefied** counts, at ASV level and at the rank set by
  `TAXA_COLLAPSE_LEVEL` (`taxa collapse`; outputs are named after the rank, e.g.
  `ancombc2_genus.qza`) `[gut_to_soil_tutorial]`. Use `DA_SAMPLE_FILTER` (an SQL `where` clause) to keep only the
  groups you compare, and `DA_REFERENCE` (e.g. `treatment::control`) to set the
  reference level.
- The visualizer is `composition da-barplot` on 2026.4+ and
  `composition ancombc2-visualizer` on 2026.1. The script switches on `Q2_EPOCH`.
- Classic ANCOM `[MH_*]` and ANCOM-BC `[moving_pictures_tutorial]` still exist. ANCOM-BC2
  adds bias correction, random effects and a pseudo-count sensitivity analysis.
- **Treat hits as hypotheses** `[gut_to_soil_tutorial]`, and confirm them with the R
  consensus (stage 11).

## Stage 10. Export and provenance (`10_export_provenance.sh`)

- Exports `feature-table.biom` / `.tsv`, `taxonomy.tsv`, `rep_seqs.fasta`, `tree.nwk` and a
  copy of the metadata. Also writes MH's combined `feature-table_w_tax.biom/.txt`. MH's
  header edit is applied to a copy, so the exported `taxonomy.tsv` is unchanged `[MH_*]`.
- `qiime tools replay-provenance --recurse` writes a re-runnable script, and
  `replay-supplement` writes a methods supplement with citations `[gut_to_soil_tutorial]`.
  (`qiime tools make-report`, mentioned in the development tutorial, is **not** in q2cli
  2026.7.)

## Stage 11. R analysis (`R/downstream_analysis.R`)

Follows `[MH_2026.3-R]`, with the defects listed in the comparison (E8–E15) fixed:

| Part | What it does | Change from MH-R |
|---|---|---|
| Settings | Reads `analysis_settings.tsv` (ranks, DA rank, overview rank, sampling depth, group, reference, seed) and prints/saves what it used | MH-R hard-coded ranks and depth |
| Import | Reads the stage-10 TSVs. Matches samples by **ID** and reports mismatches. Keeps rank names, treats empty labels with any lowercase prefix (`g__`, `sgr__`, `dvs__`) as unassigned, and **stops** if the taxonomy is deeper than the configured ranks. Prunes the SEPP reference tips from the tree. | MH-R used hard-coded column indices and got `ta1…ta7` rank names |
| QC | Library sizes and `vegan::rarecurve` with the rarefaction depth marked | — |
| Alpha | Observed, Shannon, inverse Simpson: **median of 100 rarefactions**; Wilcoxon or Kruskal–Wallis with BH correction | MH-R used one rarefaction and no test |
| Beta | Bray–Curtis (proportions), **robust Aitchison** (counts), Jaccard and unweighted UniFrac (rarefied), weighted UniFrac (proportions); `adonis2(by = "margin")` **+ `betadisper`**; PCoA, with ellipses only for categorical groups | MH-R had no dispersion test, put presence/absence metrics on unrarefied data, and drew ellipses on a continuous variable |
| Composition | Phylum- and genus-level stacked bars of the top 15 taxa + "Other" + "Unassigned" (grey) | MH-R's bars didn't sum to 100% |
| DA | Genus level, prevalence ≥10%, **same unrarefied counts** for ANCOM-BC2, ALDEx2 (Wilcoxon or KW) and MaAsLin3 (TSS/LOG + log-depth covariate). BH q ≤ 0.10. Consensus = significant in ≥2 methods **with the same direction**. ALDEx2's effect sign is aligned to the reference level, since it orders groups alphabetically. | MH-R gave MaAsLin3 rarefied input, loaded the wrong MaAsLin package, and didn't check direction |

The reference level comes from `DA_REFERENCE` (e.g. `treatment::control`). Set
`cfg$covariates` (e.g. `c("site")`) at the top of the script for adjusted models.

---

## 12. Reporting checklist

Report these in the methods section:

- [ ] QIIME 2 release, plus `replayed_commands.sh` / `reproducibility_supplement.zip` (stage 10)
- [ ] Primer sequences and the cutadapt settings
- [ ] Denoiser and its parameters (truncation and trim lengths), and the fraction of reads retained
- [ ] Rare-ASV thresholds (`04_filter_features/filter_thresholds.tsv`)
- [ ] Classifier: reference database, version, weighted or not, and the scikit-learn version
- [ ] Taxa excluded; minimum sample depth; samples removed per group
- [ ] Tree method and reference
- [ ] Even-sampling depth, number of iterations, rarefaction vs bootstrap, seed
- [ ] Statistical tests, permutations, and multiple-testing correction
- [ ] DA methods, thresholds and the consensus rule; R `sessionInfo.txt`

## 13. Validation status of this SOP

| Component | How it was checked | Result |
|---|---|---|
| Every `qiime` command and parameter | Read against plugin source at tags 2026.1.0 / 2026.4.0 / 2026.7.0 (and 2024.5.0 where relevant) | Names, defaults and outputs match |
| Bash control flow (all 11 stages) | `shellcheck` (clean) plus `tests/smoke_test.sh`: dry runs with stub `qiime`/`biom` commands on eight configurations (Illumina PE 16S DADA2 2026.7; Deblur 18S 2026.1; PacBio ITS + kmer; Illumina single-end; no primers; **18S_PR2 preset**; **preset with your own edits**; **preset with the overrides section**), with assertions on the commands issued and the settings reports. Run under **bash 5.2 and bash 4.3.48**. | All configurations pass on both bash versions. Before the empty-array fix, all five original configurations failed on bash 4.3. |
| Preset mechanism | Smoke test checks that every preset line's "shipped" value matches `config.sh`, that preset values are applied, and that values you set yourself are kept and labelled | Pass |
| R script | `tests/r_check.sh` (R 4.5, phyloseq 1.54, vegan 2.7.5, ALDEx2 1.42, maaslin3 1.2, ANCOMBC 2.14) on synthetic data with **6 known enriched genera**: 16S defaults; reference = the alphabetically second group; PR2 9 ranks at level 8; PR2 taxonomy with too few configured ranks (must stop); and a unit check of the empty-label rule | All five pass. Before the fix, the alphabetically-second reference lost all 6 true genera and the old empty-label rule left `sgr__`/`dvs__` as taxon names. |
| PR2 facts | RESCRIPt 2026.7.0 source (`get-pr2-data`, `dereplicate`); the PR2 5.1.1 taxonomy file (240,201 sequences) | 9 prefixed ranks; organelle tags `:plas :mito :apic :chro :nucl`; *Bacteriastrum* present |
| End-to-end on real data with QIIME 2 | **Not done.** QIIME 2 packages could not be installed where this was written. | Run the pipeline on a small dataset (e.g. the gut-to-soil tutorial data) before production use. |

Re-run the checks with `bash tests/smoke_test.sh <new-dir>` (no QIIME 2 needed) and
`bash tests/r_check.sh <new-dir>` (needs the R packages from §1).
