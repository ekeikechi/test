# Amplicon workflow comparison: Microbiome Helper vs. official QIIME 2 tutorials

This document compares eight sources side by side, step by step. It covers which tools,
parameters and defaults each source uses, what changed between versions, and where the
sources disagree or contain errors. The companion document
[`02-optimized-amplicon-sop.md`](02-optimized-amplicon-sop.md) turns the comparison into a
single combined SOP.

## 1. Sources and the exact revisions that were read

All sources were read from their version-controlled originals, not from rendered web pages,
so the exact revision can be stated. The amplicon-docs site links (`/en/stable/`) resolve to
the `stable` tag of each repository, which is the **QIIME 2 2026.7** release.

| Label (as given) | Page | Source read | Revision |
|---|---|---|---|
| `MH_qiime2-amplicon-2024.5` | Amplicon SOP v2 (qiime2-2024.5) | `LangilleLab/microbiome_helper.wiki` | `630a785` (2026-09-29) |
| `MH_2026.3-workflow` | Microbiome Helper 2 – Marker gene workflow | same wiki | `630a785` |
| `MH_2026.3-stats-visualization` | Microbiome Helper 2 – Amplicon basic statistics and visualisation | same wiki | `630a785` |
| `MH_2026.3-R` | Microbiome Helper 2 – Analysis workflow in R | same wiki | `630a785` |
| `moving_pictures_tutorial` | *Moving Pictures* tutorial | `qiime2/moving-pictures-tutorial` tag `stable` | `8750eba` (2026-07-21) |
| `gut_to_soil_tutorial` | Gut-to-soil axis tutorial | `caporaso-lab/gut-to-soil-tutorial` tag `stable` | `8b9601d` (2026-07-21) |
| `experienced_researchers` | QIIME 2 for experienced microbiome researchers | `qiime2/amplicon-docs` tag `stable` | `3c97730` (2026-07-21) |
| `conceptual_overview` | Conceptual overview of marker gene analysis | `qiime2/amplicon-docs` tag `stable` | `3c97730` |

> **Moving target warning.** The `main` (development) branch of `gut_to_soil_tutorial`
> (`aae3ca7`, 2026-10-08) is a substantial rewrite of the stable page (about 930 lines added
> and 330 removed). It switches the teaching classifier from Greengenes 13_8 to GTDB r202,
> switches q2-boots from bootstrapping to rarefaction (`replacement=False`), and expands
> provenance replay. The `main` branch of `moving_pictures_tutorial` adds a banner saying the
> gut-to-soil tutorial supersedes it. Where this matters, the development version is called
> out below as *g2s-dev*. Everything else describes the stable (2026.7) pages that the links
> resolve to.

### 1.1 How claims were verified

Statements about commands, parameters, defaults and version differences were checked against
the **plugin source code at the matching release tags**. Sources checked: `q2-dada2`,
`q2-deblur`, `q2-cutadapt`, `q2-vsearch`, `q2-feature-table`, `q2-feature-classifier`,
`q2-fragment-insertion`, `q2-diversity`, `q2-composition`, `q2-boots`, `q2-stats`, `q2-taxa`,
`q2-phylogeny`, `q2cli` and `rachis` at tags `2024.5.0`, `2026.1.0`, `2026.4.0` and
`2026.7.0`. Package versions in each release were checked against the
`qiime2/distributions` environment files, and classifier compatibility against
`qiime2/resources` (`53190d6`). QIIME 2 itself could not be installed in the environment
used to write this (`packages.qiime2.org` was not reachable), so no QIIME 2 command was
executed.

---

## 2. What each source is

| | MH 2024.5 | MH 2026.3 (workflow + stats + R) | Moving Pictures | Gut-to-soil (stable) | Experienced researchers / Conceptual overview |
|---|---|---|---|---|---|
| Type | Lab SOP, terse | Lab SOP and teaching workshop, split over 3 pages | Teaching tutorial | Teaching tutorial and recommended entry point | Reference explanations, no runnable commands |
| QIIME 2 version | `qiime2-amplicon-2024.5` | Activates `qiime2-amplicon-2026.1`. The text says 2025.4, and the Illumina steps are "not yet tested in 2026.1". | 2026.7 docs build | 2026.7 docs build | 2026.7 docs build |
| Example data | Your own MiSeq paired-end data | Arctic Ocean: Illumina V4-V5 **and** PacBio full-length 16S | EMP V4, **single-end**, **multiplexed** HiSeq | EMP V4, paired-end MiSeq 2×250, already demultiplexed | – |
| Platforms covered | Illumina PE | Illumina PE + PacBio CCS (HiFi) | Illumina SE | Illumina PE | All, conceptually |
| Markers covered | 16S, 18S, ITS | 16S, 18S, ITS (+12S classifier mention) | 16S | 16S | All, conceptually |
| Downstream stats | Inside QIIME 2 | QIIME 2 basics page + full R (phyloseq) page | Inside QIIME 2 | Inside QIIME 2. Group-significance tests deferred "until q2-stats is ready". | Pointers only |

---

## 3. Step-by-step comparison

Legend: **—** = step not present in that source.

### 3.1 Set-up, metadata, raw-read QC and import

| Step | MH 2024.5 | MH 2026.3 | Moving Pictures | Gut-to-soil | Explanations |
|---|---|---|---|---|---|
| Environment | conda env with FastQC + MultiQC installed into it | conda env; FastQC/MultiQC installed separately | amplicon distribution | any install; `qiime info` check | – |
| Metadata | Format in spreadsheet; Keemei validation; `$METADATA` variable | Same. Tutorial-specific `cut`/`sed` fix-ups. Metadata **not used** on the workflow page. | `metadata tabulate` is the first command | `metadata tabulate` is the first command; merges DADA2 stats into metadata | TSV; first column `sample-id`/`id`; never include PII (MP) |
| Raw-read QC | FastQC + MultiQC (per-base Q > 30 rule of thumb) | Same, plus `demux summarize` on the **raw** import | `demux summarize` | `demux summarize` | – |
| Import | `CasavaOneEightSingleLanePerSampleDirFmt` (strict file naming) | Casava **or** manifest (`PairedEndFastqManifestPhred33V2` / `SingleEndFastqManifestPhred33V2`); lab `make_manifest.py` helper | `EMPSingleEndSequences` + `demux emp-single` | Pre-imported artifact; recommends a manifest / asking the sequencing centre | "Importing with a manifest file will always work" |
| Demultiplexing | Assumes demultiplexed | Assumes demultiplexed | `demux emp-single` (Golay barcodes) | Already demultiplexed | `emp-paired`, `cutadapt demux-*`; dual-index PE must be done outside QIIME 2 |

### 3.2 Primer removal and denoising

| Step | MH 2024.5 | MH 2026.3 | Moving Pictures | Gut-to-soil | Explanations |
|---|---|---|---|---|---|
| Primer removal | `cutadapt trim-paired --p-discard-untrimmed --p-no-indels`. Default example is **16S V6–V8** (bacteria). Alternatives: V4–V5, archaeal V6–V8, 18S V4, ITS2. | Same command. Default example changed to **16S V4–V5** ("universal"). Link to the IMR primer list. PacBio: primers removed **inside** `dada2 denoise-ccs`. | — (EMP reads lack primers) | — on stable. g2s-dev: not needed for EMP; use DADA2 trim-left if at a fixed position, cutadapt otherwise. | Use cutadapt, or DADA2 `--p-trim-left*` |
| Denoiser (Illumina) | **Deblur** (DADA2 on a separate quick-reference page) | **Deblur** | DADA2 `denoise-single` **or** Deblur | **DADA2 `denoise-paired`** | DADA2 or Deblur; clustering discouraged |
| Read merging | `vsearch merge-pairs` (defaults; "more stringent than PEAR") | Same | n/a (single-end) | Done inside DADA2 | "Do not merge" for DADA2; merge first for Deblur or OTUs |
| Pre-denoise quality filter | `quality-filter q-score` (defaults) | Same | `quality-filter q-score` before Deblur only | none | Needed for Deblur, unnecessary for DADA2 |
| Length / truncation | `deblur --p-trim-length -1` (raise it only if Deblur errors) | Same. New guidance on choosing trim length from the length-percentile table (e.g. 366 nt for V4–V5). | DADA2 `trunc-len 120` / Deblur `trim-length 120` | DADA2 trim-left 0/0, trunc-len 250/250 | Optional; reads shorter than the truncation length are dropped |
| Deblur options | `denoise-16S --p-sample-stats --p-min-reads 1`. 18S/ITS: `denoise-other` with PR2 / UNITE positive-filter references. | Same | `denoise-16S --p-sample-stats` (default `min-reads` 10) | n/a | `denoise-16S` positive-filters against GG 13_8 85% OTUs (≥60% identity) |
| Denoiser (PacBio) | — | `dada2 denoise-ccs` with `--p-front/--p-adapter`, `--p-min-len 1200 --p-max-len 1800` (16S), `--p-max-ee 18`, `--p-pooling-method pseudo`, `--p-n-reads-learn 100000`. 18S and ITS variants. | — | — | `denoise-ccs` for PacBio CCS |
| Denoising QC | `feature-table summarize`. Check that depth did not collapse after Deblur. | Same; outputs copied to a common `denoised_output/` | `metadata tabulate` on DADA2 stats; `deblur visualize-stats` | `metadata tabulate` on DADA2 stats, merged with sample metadata | – |
| OTU clustering | — | — | — | — (g2s-dev: explained as discouraged) | Possible with q2-vsearch; **recommended against** |

### 3.3 Feature table summaries, filtering and taxonomy

| Step | MH 2024.5 | MH 2026.3 | Moving Pictures | Gut-to-soil | Explanations |
|---|---|---|---|---|---|
| Table summary | `feature-table summarize --o-visualization` | `feature-table summarize --o-summary --o-sample-frequencies --o-feature-frequencies` (new pipeline interface) | `summarize` + `tabulate-seqs` | `summarize` (with metadata) + `tabulate-seqs` with feature frequencies as metadata | `summarize` + `tabulate-seqs` together |
| Abundance / prevalence filter | `filter-features --p-min-frequency` = round(0.001 × mean sample depth) ("MiSeq bleed-through"), `--p-min-samples 1` | Same, with the note that the 0.1% rule may be outdated | none | `filter-features --p-min-samples 2` (drop single-sample ASVs) | Filter low-abundance features **before** classifying to save time |
| Order of filter vs. taxonomy | Classify first, then filter | Same | Classify unfiltered | **Filter first**, then classify only the retained sequences | Same advice as g2s |
| Classifier | sklearn NB, SILVA 138 99% full-length (2022.11-era file path); custom UNITE v10 ITS classifiers | sklearn NB. Downloads for **GG2 2024.09**, **SILVA 138**, **GTDB r220** (all sklearn 1.4.2). Lab mirror. HOMD build example. Region-specific classifiers deliberately **not** used. | Trains an NB classifier on GG 13_8 85% OTUs, labelled "suboptimal" for teaching | Same as MP on stable; g2s-dev trains on GTDB r202 via RESCRIPt `get-gtdb-data` with reduced `n_features` | NB "can slightly outperform" alignment methods; **environment-weighted classifiers recommended**; 8–16+ GB RAM |
| `classify-sklearn` options | `--p-n-jobs` | adds `--p-read-orientation same` | defaults | defaults (confidence 0.7) | – |
| Taxonomy QC | `tabulate-seqs`; BLAST ~5 ASVs from different phyla | Same, plus `grep` of `taxonomy.tsv` | BLAST spot-check | `tabulate-seqs` with the taxonomy collection; BLAST spot-check; trust the curated classifier over BLAST | `tabulate-seqs` to compare classifiers |
| Taxonomy-based filtering | `taxa filter-table --p-exclude mitochondria,chloroplast`. Optional `--p-include p__`, which is **not** recommended by default. | Same; skip for non-16S markers | — | — | Mentions contaminant/host filtering |
| Low-depth samples | Look at `observed_features` rarefaction curves, then `filter-samples --p-min-frequency` (typically 1000–4000), or keep all | Same, plus "try several cut-offs" | Handled implicitly by the rarefaction depth | Remove obviously failed samples with `filter-samples`; check per-group retention | – |

### 3.4 Phylogeny and diversity

| Step | MH 2024.5 | MH 2026.3 | Moving Pictures | Gut-to-soil | Explanations |
|---|---|---|---|---|---|
| 16S tree | **SEPP** insertion, reference `sepp-refs-gg-13-8` | **SEPP** insertion, reference **GG2 2022.10 backbone** (or SILVA 128 / GG 13_8) | `align-to-tree-mafft-fasttree` (de novo) | **No tree**: kmer-based diversity instead | Reference-based trees are "generally higher quality"; de novo trees "often fairly low-quality" |
| 18S/ITS tree | de novo: `mafft` → `mask` → `fasttree` → `midpoint-root` (4 commands) | Same 4 commands | n/a | n/a | – |
| Rarefaction curves | `alpha-rarefaction` twice: with metadata (grouped) and without (per sample); `max-depth` = maximum sample depth | Stats page: once, without metadata | `max-depth 4000` (≈ median); with metadata | `max_depth 260` (≈ median); with metadata | – |
| Diversity engine | `core-metrics-phylogenetic` (single rarefaction) | Same (stats page) | `core-metrics-phylogenetic`, depth 1103 | **q2-boots `kmer-diversity`**: depth 96, `n=10`, **bootstrapping** (`replacement=True`), median/medoid averaging. g2s-dev: rarefaction, `n=100` recommended. | q2-boots "results are more robust" |
| Alpha stats | `alpha-group-significance` (Shannon) | Same, for each metric | `alpha-group-significance` (Faith PD, evenness); mentions `alpha-correlation` | — (deferred to q2-stats) | – |
| Beta stats | — | `beta-group-significance` (Bray–Curtis PERMANOVA) | `beta-group-significance --p-pairwise` (unweighted UniFrac); mentions `mantel`, `bioenv` | — | – |
| Ordination | Emperor (from core-metrics) | Emperor | Emperor with `--p-custom-axes days-since-experiment-start` | q2-vizard 2-D scatter plot | – |

### 3.5 Composition, differential abundance and export

| Step | MH 2024.5 | MH 2026.3 | Moving Pictures | Gut-to-soil | Explanations |
|---|---|---|---|---|---|
| Barplots | `taxa barplot`, plus `feature-table group --p-mode sum` per category | `taxa barplot` | `taxa barplot` | `taxa barplot` (always on the ASV table, never the kmer table) | `taxa barplot`, `feature-table heatmap` |
| Differential abundance (QIIME 2) | **ANCOM** (`add-pseudocount 1` → `composition ancom`) | Same (stats page) | **ANCOM-BC**: gut samples only, `formula subject`; ASV and L6 (genus) | **ANCOM-BC2**: 3 sample types, reference level = HEC; ASV and genus (`taxa collapse --p-level 6`) | ANCOM-BC; "consult a statistician" |
| Differential abundance (R) | Points to ALDEx2 / corncob | **ALDEx2 + MaAsLin3 + radEmu**, consensus of taxa found by ≥2 (or all 3) tests | — | — | – |
| Export | FASTA; BIOM + taxonomy → TSV (`biom add-metadata` / `biom convert`); Newick tree | Same | — | **Provenance replay** to regenerate a script | `qiime tools export` |
| Further analyses | PICRUSt2, STAMP, Random Forest, corncob | Phylo-RPCA (gemelli) | — | sample-classifier, longitudinal (as exercises) | longitudinal, sample-classifier, q2-kmerizer, q2-boots, q2-micom |

### 3.6 R workflow (`MH_2026.3-R`), the only source that leaves QIIME 2

| Aspect | What it does | Comment |
|---|---|---|
| Import | `read.csv` on `feature-table_w_tax.txt`; taxonomy split into 7 ranks; metadata columns picked by index; `read_tree`; `phyloseq()` | Fragile (hard-coded indices); see errata E10–E11 |
| Filtering | `prune_taxa(taxa_sums ≥ 100)`, `prune_samples(sample_sums ≥ 1000)` | A **second**, different filtering rule on top of the QIIME 2 0.1% rule |
| Normalisation | Rarefy once to the minimum depth (seed 711) for alpha diversity; relative abundance (%) for beta diversity and plots | Unweighted metrics on non-rarefied data are depth-sensitive |
| Alpha | `plot_richness` (Observed, Shannon, Fisher, Simpson) | No test reported |
| Beta | Bray–Curtis and weighted UniFrac; `adonis2(..., by = "margin")`; PCoA + `stat_ellipse` | No dispersion check; ellipses drawn on a continuous variable |
| Composition | Top-20 ASV / phylum / genus barplots; `plot_heatmap`; `plot_tree` | No "Other" bar, so stacks do not reach 100% |
| DA | Genus level: ALDEx2 (`aldex.kw`), MaAsLin3 (`TSS` + `LOG`, `~ group + sample_sums`), radEmu (`emuFit`); keep taxa found by all three, plot | Good multi-method idea; implementation issues in E8, E9 and E12 |
| Bonus | Merges Illumina and PacBio results at species/genus level | See E9 |

---

## 4. What changed: MH 2024.5 → MH 2026.3

1. **Split into pages.** Upstream processing, QIIME 2 statistics and R analysis are now
   separate; metadata is no longer used on the upstream page.
2. **PacBio HiFi path added**, using `dada2 denoise-ccs` with primer removal, orientation
   fixing and length gating inside DADA2, for 16S (bacteria, archaea), 18S and ITS.
3. **Manifest import** (and the lab's `make_manifest.py`) added alongside Casava import, plus
   a raw-read `demux summarize`.
4. **Default primer example** changed from 16S V6–V8 to 16S V4–V5.
5. **Deblur trim-length guidance**: pick a length from the length-percentile table instead of
   leaving `-1` by default.
6. **`feature-table summarize`** updated to the 2026 pipeline interface (three outputs).
7. **Unified file names** (`denoised_output/table.qza`, `rep_seqs_final.qza`) so the Deblur and
   DADA2 branches feed one downstream path.
8. **Taxonomy**: a choice of seven databases. Direct links to the QIIME 2 sklearn-1.4.2 builds of
   GG2 2024.09, SILVA 138 and GTDB r220. A HOMD custom-classifier example.
   `--p-read-orientation same`.
9. **SEPP reference** changed from Greengenes 13_8 to the Greengenes2 2022.10 backbone.
10. **Downstream statistics** gained `beta-group-significance` (PERMANOVA); ANCOM is kept.
11. **New R workflow** with phyloseq/vegan plots and multi-tool differential abundance
    (ALDEx2, MaAsLin3, radEmu).
12. **Removed from the upstream page:** grouped barplots (`feature-table group`) and the
    with-metadata/without-metadata pair of rarefaction plots.

## 5. What changed in QIIME 2 itself (verified against source), and why it matters here

| Change | 2024.5 | 2026.1 | 2026.4 / 2026.7 | Consequence for these workflows |
|---|---|---|---|---|
| Distribution name / env | `qiime2-amplicon-2024.5` | `qiime2-amplicon-2026.1` | **`qiime2`** distribution (env file `rachis-qiime2-linux-64-conda.yml`; Library quickstart env name `rachis-qiime2-<epoch>`) | `conda activate qiime2-amplicon-…` lines in MH need updating for 2026.4+ |
| scikit-learn | 1.4.2 | 1.4.2 | **1.7.1** | `classify-sklearn` raises `ValueError` when the classifier's sklearn version differs. MH 2026.3 classifier links (sklearn 1.4.2) work on 2024.5–2026.1 and **fail on 2026.4+**. For 2026.4+ the QIIME 2 resources page points to SILVA (arb-silva.de) and GTDB (Zenodo `10.5281/zenodo.19446101`) builds. |
| `feature-table summarize` | Visualizer (`--o-visualization`); `summarize-plus` pipeline exists | **Pipeline**: `--o-summary`, `--o-sample-frequencies`, `--o-feature-frequencies` | same as 2026.1 | MH 2024.5 summarize commands fail on 2026.x |
| `dada2 denoise-*` outputs | table, rep-seqs, stats | adds `--o-base-transition-stats` | same; `denoise-paired` also gains `--p-retain-unmerged` | Every DADA2 call needs the extra output (or `--output-dir`) |
| q2-boots, q2-kmerizer, q2-stats, q2-vizard | not in the distribution | **included** | included | `conceptual_overview` ("isn't yet part of the amplicon distribution") is outdated |
| ANCOM-BC2 | not available (ANCOM-BC 2.4.0 only) | `composition ancombc2`; plot with `composition ancombc2-visualizer` | `ancombc2`; plot with **`composition da-barplot`** (accepts both ANCOM-BC and ANCOM-BC2 output) | The visualizer name depends on the release |
| Classic ANCOM | `composition ancom` | still present | still present | MH stats page still runs, but QIIME 2 tutorials have moved to ANCOM-BC / ANCOM-BC2 |
| q2-stats `alpha-group-significance` | – | – | present but **deprecated** | Use `diversity alpha-group-significance` for now |
| `diversity adonis` | present | present | present; runs `adonis2(by = 'terms')` (sequential) | Not directly comparable with MH-R's `by = "margin"` |
| `qiime tools make-report` | – | – | **not present** in q2cli 2026.7 (only mentioned by g2s-dev) | Don't rely on it with 2026.7 |

---

## 6. Errata and inconsistencies found in the sources

Severity: **High** = a command will fail or the result is wrong; **Medium** = a likely silent
problem or misleading result; **Low** = clarity or typo.

| # | Source | Issue | Evidence | Severity |
|---|---|---|---|---|
| E1 | MH 2026.3 | Classifier downloads are sklearn-1.4.2 builds. They fail on QIIME 2 2026.4+ (sklearn 1.7.1). | `q2_feature_classifier/_taxonomic_classifier.py` raises `ValueError` on version mismatch; distribution env files | High (on 2026.4+) |
| E2 | MH 2024.5 | Classifier path `…/qiime2-2022.11_classifiers/silva-138-99-nb-classifier.qza` suggests a 2022.11-era build (sklearn 0.24.1 per the resources page). That fails on 2024.5 (sklearn 1.4.2) unless it was retrained locally. | resources index "QIIME 2 2021.4–2024.2 … Sklearn 0.24.1" | Medium (likely) |
| E3 | MH 2024.5 | `feature-table summarize --o-visualization` does not exist on 2026.x | q2-feature-table 2026.1.0 registers `summarize` as a pipeline | High (on 2026.x) |
| E4 | MH 2026.3 | Version statements disagree: text says 2025.4, the command activates 2026.1, and the Illumina steps are "not yet tested" | page text | Low |
| E5 | MH 2026.3 | PacBio text says the length window is 1300–1800 nt; the command uses `--p-min-len 1200` | page text vs. command | Low |
| E6 | MH both | "`--p-min-reads 1` will retain singletons" is only partly true. It relaxes the across-sample threshold (default 10), but `--p-min-size` (default 2) still drops per-sample singletons. | q2-deblur parameter descriptions | Low |
| E7 | MH both | SEPP output is used directly. Fragments SEPP cannot place are absent from the tree, and `fragment-insertion filter-features` (which exists for exactly this) is never run, so phylogenetic metrics can fail or drop data. | q2-fragment-insertion 2026.7.0 `filter_features` description | Medium |
| E8 | MH R | Loads `library(Maaslin2)` but calls `maaslin3()`; never loads `radEmu` (`emuFit`) | page code | High |
| E9 | MH R | `merge_phyloseq(physeq_illum, physeq_pacbio)` merges the **ASV-level** objects, not the species-glommed `physeq_species_*`. The glom and rename steps are dead code and the merged object contradicts the stated intent. | page code | High (logic) |
| E10 | MH R | `tax_table(<data.frame>)` drops column names, so ranks become `ta1…ta7` (hence `"ta2"`, `"ta6"`). `as.matrix()` keeps `Phylum`, `Genus`, … | phyloseq `taxonomyTable-class.R` | Low |
| E11 | MH R | Hard-coded column indices (`[,2:19]`, `c(3,4,5,6,7)`, `c(2,4,5)`) and hand-typed group vectors in sample order. Silent mislabelling on any other dataset. | page code | Medium |
| E12 | MH R | MaAsLin3 is run on the **rarefied** genus table, while ALDEx2 and radEmu use unrarefied counts. The consensus therefore mixes inputs. | page code | Medium |
| E13 | MH R | `stat_ellipse` on a continuous colour variable (longitude) | page code | Low |
| E14 | MH R / QIIME 2 | PERMANOVA without a dispersion (`betadisper`) check. MH-R uses `by = "margin"` while QIIME 2 `adonis` uses `by = "terms"`, so p-values and R² can differ for correlated predictors. | `q2_diversity/_beta/adonis_assets/run_adonis.R` | Medium |
| E15 | MH R | Top-N barplots without an "Other" category, so bars don't sum to 100% | page code | Low |
| E16 | MH 2026.3 | DADA2 "copy" section says "we'll just copy the Deblur output" | page text | Low |
| E17 | MH 2024.5 | Table-of-contents links for 3.1–3.3 point to the 2020.8 page | page text | Low |
| E18 | conceptual_overview | Says q2-boots is not in the amplicon distribution. It has been since 2026.1. | distribution env files | Low |
| E19 | experienced_researchers | Lists `denoise-css` (should be `denoise-ccs`) and links it to the `denoise-pyro` anchor | page text | Low |
| E20 | gut_to_soil (stable) | `trim-left-r` / `trunc-len-r` described as acting on "each forward read" | page text | Low |

---

## 7. Where the sources agree, where they diverge, and which choice is better

**Agreement (strong defaults):**

- Denoise into ASVs. Don't cluster into OTUs (`conceptual_overview`, `experienced_researchers`,
  g2s-dev).
- Use an sklearn Naive Bayes classifier trained on a **full-length** reference (MH 2026.3,
  resources page), version-matched to your QIIME 2. Spot-check with BLAST.
- Remove mitochondria and chloroplast for 16S. Don't blanket-remove unclassified features
  (MH both).
- Pick an even-sampling depth that balances samples lost against sequences kept, guided by
  alpha-rarefaction curves (all).
- Treat differential-abundance hits as hypotheses, and use compositional methods (all; MH-R
  goes furthest with a multi-tool consensus).

**Divergence, and the recommended resolution** (adopted in the SOP):

| Topic | Options in the sources | Recommended | Why |
|---|---|---|---|
| Illumina denoiser | Deblur (MH) vs DADA2 (g2s, MP option) | **DADA2 `denoise-paired`**, with Deblur kept as an option for 16S | DADA2 needs no pre-merge or q-score filter (`experienced_researchers`). It avoids the read loss from stringent `vsearch merge-pairs` that MH itself warns about. It works for 18S/ITS without custom positive-filter references. The same tool is used for PacBio. |
| Filter order | classify → filter (MH) vs filter → classify (g2s, conceptual) | **Abundance/prevalence filter → classify → taxonomy filter** | Faster and lower memory; taxonomy filtering still happens after classification |
| Rare-feature filter | 0.1% of mean depth (MH) vs ≥2 samples (g2s) vs ≥100 reads (MH-R) | **One** configurable rule in QIIME 2 (default: ≥2 samples **and** ≥0.1% of mean depth); no second filter in R | Double filtering with different rules (MH QIIME + MH R) is hard to report and inconsistent |
| Tree | SEPP (MH) vs de novo (MP) vs none/kmer (g2s) | **SEPP + `filter-features`** for 16S; `align-to-tree-mafft-fasttree` for 18S; kmer-based or non-phylogenetic metrics for ITS | Reference placement is higher quality for short 16S reads (`experienced_researchers`). ITS is not alignable across fungi. |
| Diversity | single rarefaction (MH, MP) vs repeated resampling (g2s) | **q2-boots `core-metrics`** (rarefaction, n = 100; phylogeny when available) | Uses all the data and is less sensitive to one random draw (g2s, conceptual) |
| Beta-diversity stats | none (g2s), PERMANOVA (MH-stats, MP), adonis2 (MH-R) | PERMANOVA pairwise + multi-factor `adonis`, **plus a dispersion check** | PERMANOVA is sensitive to unequal dispersion |
| Differential abundance | ANCOM (MH), ANCOM-BC (MP), ANCOM-BC2 (g2s), ALDEx2 + MaAsLin3 + radEmu (MH-R) | **ANCOM-BC2** in QIIME 2, **plus an R consensus** (ANCOM-BC2 + ALDEx2 + MaAsLin3, all on the same unrarefied counts, ≥2 of 3) | Combines the current QIIME 2 recommendation with MH-R's multi-method practice, without MH-R's input inconsistency |
| Reproducibility | "keep your commands" (MH) vs provenance replay (g2s) | **Both**: scripted pipeline + `replay-provenance` + `replay-supplement` | – |
