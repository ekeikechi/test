# Amplicon workflow: comparison and combined SOP

This repository compares the Microbiome Helper amplicon SOPs with the official QIIME 2
amplicon documentation, then combines them into one optimized, runnable workflow
(QIIME 2 CLI + R).

| File | What it is |
|---|---|
| [`docs/01-workflow-comparison.md`](docs/01-workflow-comparison.md) | Step-by-step comparison of all eight sources: tools, parameters and defaults; what changed between MH 2024.5 and MH 2026.3 and across QIIME 2 releases; 20 errata with evidence; recommended resolutions |
| [`docs/02-optimized-amplicon-sop.md`](docs/02-optimized-amplicon-sop.md) | The combined SOP: decisions and rationale, installation, eleven stages with checkpoints, reporting checklist, validation status. Every step cites its source label. |
| [`pipeline/`](pipeline) | `config.sh` (all settings and decisions, plus an `18S_PR2` preset), `lib.sh`, stage scripts `01_…10_`, optional `05a_build_pr2_classifier.sh`, and `auto_truncation.py` (automatic DADA2 truncation), `denoise_check.py` (post-DADA2 sample flags) and `choose_depth.py` (minimum, automatic and sensitivity rarefaction depths) (QIIME 2 2026.7; 2026.1 supported; bash 3.2+) |
| [`R/downstream_analysis.R`](R/downstream_analysis.R) | Stage 11: phyloseq/vegan diversity, composition plots, and an ANCOM-BC2 + ALDEx2 + MaAsLin3 differential-abundance consensus, for independent samples |
| [`R/design_analysis.R`](R/design_analysis.R), [`R/designs/`](R/designs) | Stage 11b, nested and split-plot designs: linear mixed models for alpha diversity, hierarchical PERMANOVA with restricted permutations, ANCOM-BC2 + MaAsLin3 with plot and subplot random effects (ALDEx2 supporting); one design file per analysis, including the two Grande Anse manuscripts |
| [`R/pr2_traits.R`](R/pr2_traits.R), [`R/pr2_trait_overrides.tsv`](R/pr2_trait_overrides.tsv), [`R/adl2019_trophic.tsv`](R/adl2019_trophic.tsv) | Stage 11a (18S + PR2): every PR2 trait for each ASV (ecological function, mixoplankton type, HAB, WoRMS habitat, ids, source evidence), dinoflagellate photosynthesis and diets for all groups from Adl et al. 2019 Appendix S2, a photosynthesis category, cited overrides, a habitat check, and exports of the photosynthetic set (main and sensitivity) |
| [`tests/`](tests) | `smoke_test.sh` (12 pipeline configurations against stub `qiime`/`biom`, with assertions on commands, presets and settings reports; the stub rejects commands whose interface differs in the configured QIIME 2 release), `r_check.sh` (R script on synthetic 16S and PR2 data with known answers, plus negative checks), `auto_truncation_test.sh` (automatic DADA2 truncation on synthetic quality profiles), `traits_check.sh` (PR2 trait annotation on a synthetic PR2 database with known answers), `qc_scripts_test.sh` (post-DADA2 flags and depth rules on hand-computed cases), `design_check.sh` (split-plot statistics on a synthetic 216-sample experiment with known effects in each stratum) |

## Sources (labels used throughout)

| Label | URL |
|---|---|
| `MH_qiime2-amplicon-2024.5` | https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2%E2%80%902024.5) |
| `MH_2026.3-workflow` | https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-Helper-2-Marker-gene-workflow |
| `MH_2026.3-stats-visualization` | https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-Helper-2-Amplicon-basic-statistics-and-visualisation |
| `MH_2026.3-R` | https://github.com/LangilleLab/microbiome_helper/wiki/Microbiome-Helper-2-Analysis-workflow-in-R |
| `moving_pictures_tutorial` | https://amplicon-docs.qiime2.org/en/stable/tutorials/moving-pictures/ |
| `gut_to_soil_tutorial` | https://amplicon-docs.qiime2.org/en/stable/tutorials/gut-to-soil/ |
| `experienced_researchers` | https://amplicon-docs.qiime2.org/en/stable/explanations/experienced-researchers/ |
| `conceptual_overview` | https://amplicon-docs.qiime2.org/en/stable/explanations/conceptual-overview/ |

Exact revisions read, and how claims were verified against plugin source code, are listed in
section 1 of the comparison.

## Quick start

```bash
conda activate rachis-qiime2-2026.7          # see SOP §1 for installation
$EDITOR pipeline/config.sh                   # paths, metadata columns; PRESET="18S_PR2" for 18S V4 + PR2
bash pipeline/01_import_qc.sh                # each stage stops at a checkpoint when it needs a decision
# ... 02 to 10 ...
Rscript R/downstream_analysis.R results/10_export results/11_R   # reads ranks, depth, group from stage 10
Rscript R/pr2_traits.R results/10_export results/11_traits           # 18S + PR2: traits, photosynthesis category
Rscript R/design_analysis.R R/designs/ms2_community_all_depths.tsv   # split-plot statistics (SOP §11b)
```

## Headline findings from the comparison

- The classifier links in **MH 2026.3 are scikit-learn 1.4.2 builds**. They work on QIIME 2
  2024.5–2026.1 but **fail on 2026.4+** (scikit-learn 1.7.1), where QIIME 2 now points to
  SILVA- and GTDB-hosted classifiers.
- **`feature-table summarize`** changed interface in 2026.1, so the MH 2024.5 commands fail on
  current releases. MH 2026.3 is already updated.
- MH's SEPP step leaves out **`fragment-insertion filter-features`**, which is needed when
  SEPP rejects any ASV.
- The MH R workflow has several bugs: wrong MaAsLin package loaded, an ASV-level merge where a
  species-level one was intended, MaAsLin3 given rarefied input while the other tools get
  raw counts, and hard-coded column indices. All are fixed in `R/downstream_analysis.R`.
- QIIME 2's own recommendations have moved to DADA2 without pre-merging, repeated-rarefaction
  diversity (q2-boots), and ANCOM-BC2. The SOP adopts all three and keeps MH's strengths:
  primer screening, contaminant filtering, SEPP, PacBio HiFi parameters, and multi-tool
  differential abundance.

## Validation status

The bash pipeline was checked with `shellcheck` and dry-run through twelve configurations
against command stubs, on bash 5.2 and 4.3. Every `qiime` command and parameter was read against the plugin source
at the matching release tags. The R script was run end to end on synthetic data, where it
recovered exactly the known differentially abundant genera. The split-plot statistics were
run on a synthetic experiment with the Grande Anse layout. Each known effect was found in
its own stratum, with the correct degrees of freedom and no false-positive DA calls. The
first run on real data
(QIIME 2 2026.7, 18S V4, 216 samples) completed stage 1. Stage 2 then stopped on an interface
change missed in the source review: `cutadapt trim-paired` in 2026.7 requires `--o-stats`. A
re-check of every plugin the pipeline calls, comparing 2026.4 with 2026.7, found one more:
the ANCOM-BC2 visualizer, used in stage 9. Both are fixed, and the test stub now rejects both
mistakes. Stage 3 then failed in `feature-table summarize` on metadata with sample IDs only;
this is fixed. Stages 1–4 have completed on real data. Stages 5–11 have not yet run on real
data, so do a small test run before production use.
