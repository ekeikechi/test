# Amplicon workflow: comparison and combined SOP

This repository compares the Microbiome Helper amplicon SOPs with the official QIIME 2
amplicon documentation, then combines them into one optimized, runnable workflow
(QIIME 2 CLI + R).

| File | What it is |
|---|---|
| [`docs/01-workflow-comparison.md`](docs/01-workflow-comparison.md) | Step-by-step comparison of all eight sources: tools, parameters and defaults; what changed between MH 2024.5 and MH 2026.3 and across QIIME 2 releases; 20 errata with evidence; recommended resolutions |
| [`docs/02-optimized-amplicon-sop.md`](docs/02-optimized-amplicon-sop.md) | The combined SOP: decisions and rationale, installation, eleven stages with checkpoints, reporting checklist, validation status. Every step cites its source label. |
| [`pipeline/`](pipeline) | `config.sh` (all settings and decisions), `lib.sh`, and stage scripts `01_…10_` (QIIME 2 2026.7; 2026.1 supported) |
| [`R/downstream_analysis.R`](R/downstream_analysis.R) | Stage 11: phyloseq/vegan diversity, composition plots, and an ANCOM-BC2 + ALDEx2 + MaAsLin3 differential-abundance consensus |
| [`tests/`](tests) | `smoke_test.sh` (all pipeline branches against stub `qiime`/`biom`), `r_check.sh` (R script on synthetic data with known answers) |

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
$EDITOR pipeline/config.sh                   # paths, platform, marker, primers, metadata columns
bash pipeline/01_import_qc.sh                # each stage stops at a checkpoint when it needs a decision
# ... 02 to 10 ...
Rscript R/downstream_analysis.R results/10_export results/11_R <group_column>
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

The bash pipeline was checked with `shellcheck` and dry-run through five configurations
against command stubs. Every `qiime` command and parameter was read against the plugin source
at the matching release tags. The R script was run end to end on synthetic data, where it
recovered exactly the known differentially abundant genera. **The pipeline has not been run
on real data with QIIME 2 installed** (QIIME 2 packages were unreachable from the environment
it was written in). Do a small test run before production use.
