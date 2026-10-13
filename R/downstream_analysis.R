#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Downstream statistics and figures for the combined amplicon SOP (stage 11).
#
# Inputs : the plain-text exports written by pipeline/10_export_provenance.sh
#          (feature-table.tsv, taxonomy.tsv, metadata.tsv, optional tree.nwk) and
#          analysis_settings.tsv, which carries TAXONOMY_RANKS, TAXA_COLLAPSE_LEVEL,
#          OVERVIEW_RANK_LEVEL, SAMPLING_DEPTH, GROUP_COL, DA_REFERENCE and RANDOM_SEED
#          from pipeline/config.sh.
# Usage  : Rscript R/downstream_analysis.R [export_dir] [out_dir] [group_column]
#          (group_column overrides GROUP_COL from the settings file)
#
# Structure follows [MH_2026.3-R] (phyloseq import, QC, alpha/beta diversity,
# composition plots, multi-tool differential abundance), with the fixes listed in
# docs/01-workflow-comparison.md section 6 (E8-E15):
#   * no hard-coded column indices; samples are matched by ID and mismatches reported
#   * taxonomy keeps real rank names (as.matrix before tax_table)
#   * alpha diversity averaged over repeated rarefactions [gut_to_soil_tutorial: q2-boots]
#   * PERMANOVA (marginal terms) reported together with a dispersion test
#   * every DA tool gets the same unrarefied counts, BH-adjusted q-values, and the
#     consensus also requires agreement on effect direction
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  export_dir     = if (length(args) >= 1) args[[1]] else "results/10_export",
  out_dir        = if (length(args) >= 2) args[[2]] else "results/11_R",
  group          = if (length(args) >= 3) args[[3]] else NULL,  # else GROUP_COL from settings
  group_ref      = NULL,        # from DA_REFERENCE ("column::level"); NULL = first level alphabetically
  covariates     = character(0),# extra adjustment terms, e.g. c("site")
  # Defaults below apply only when analysis_settings.tsv (written by stage 10) is absent.
  ranks          = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  collapse_level = 6,           # TAXA_COLLAPSE_LEVEL: rank number used for DA
  overview_level = 2,           # OVERVIEW_RANK_LEVEL: rank number for overview barplots
  rarefy_depth   = NULL,        # SAMPLING_DEPTH; NULL = smallest library size
  rarefy_iter    = 100,         # iterations averaged for alpha diversity
  seed           = 711,         # RANDOM_SEED [MH_2026.3-R]
  permutations   = 999,
  top_n          = 15,          # taxa shown in barplots; the rest are pooled as "Other"
  da_prevalence  = 0.10,        # keep taxa present in >= 10% of samples for DA
  da_q           = 0.10,        # FDR threshold used for every DA method [MH_2026.3-R]
  da_min_methods = 2            # consensus: significant in >= 2 methods [MH_2026.3-R]
)
say <- function(...) message(sprintf(...))

# Shared settings from pipeline/config.sh, handed over by stage 10. This makes
# config.sh the single source for ranks, DA rank, sampling depth, group, reference
# and seed, so QIIME 2 and R cannot silently disagree.
src <- setNames(rep("script default", 7),
                c("ranks", "collapse_level", "overview_level", "rarefy_depth", "seed",
                  "group", "group_ref"))
if (!is.null(cfg$group)) src[["group"]] <- "command line"
settings_path <- file.path(cfg$export_dir, "analysis_settings.tsv")
if (file.exists(settings_path)) {
  st <- read.delim(settings_path, colClasses = "character", quote = "",
                   na.strings = character(0))
  st <- setNames(st$value, st$setting)
  get <- function(k) if (k %in% names(st) && nzchar(st[[k]])) st[[k]] else NULL
  from_file <- function(field, value) { src[[field]] <<- "analysis_settings.tsv"; value }
  if (!is.null(v <- get("taxonomy_ranks")))
    cfg$ranks <- from_file("ranks", strsplit(v, "\\s+")[[1]])
  if (!is.null(v <- get("collapse_level"))) cfg$collapse_level <- from_file("collapse_level", as.integer(v))
  if (!is.null(v <- get("overview_level"))) cfg$overview_level <- from_file("overview_level", as.integer(v))
  if (!is.null(v <- get("sampling_depth"))) cfg$rarefy_depth <- from_file("rarefy_depth", as.integer(v))
  if (!is.null(v <- get("random_seed")))    cfg$seed <- from_file("seed", as.integer(v))
  if (is.null(cfg$group) && !is.null(v <- get("group_col"))) cfg$group <- from_file("group", v)
  if (!is.null(v <- get("da_reference")) && !is.null(cfg$group)) {
    # DA_REFERENCE may hold several "column::level" items separated by spaces
    hit <- grep(paste0("^", cfg$group, "::"), strsplit(v, "\\s+")[[1]], value = TRUE)
    if (length(hit)) cfg$group_ref <- from_file("group_ref", sub("^[^:]*::", "", hit[1]))
  }
} else {
  say("No analysis_settings.tsv in %s; using the defaults in this script.", cfg$export_dir)
}
if (is.null(cfg$group)) cfg$group <- "treatment"
for (lv_name in c("collapse_level", "overview_level"))
  if (cfg[[lv_name]] < 1 || cfg[[lv_name]] > length(cfg$ranks))
    stop(sprintf("%s = %d but only %d ranks are defined: %s", lv_name, cfg[[lv_name]],
                 length(cfg$ranks), paste(cfg$ranks, collapse = " ")))
cfg$da_rank       <- cfg$ranks[cfg$collapse_level]
cfg$overview_rank <- cfg$ranks[cfg$overview_level]

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(...) file.path(cfg$out_dir, ...)

# Report the settings actually used, with their source.
used <- data.frame(
  setting = c("ranks", "collapse_level (DA rank)", "overview_level", "rarefy_depth", "seed",
              "group", "group_ref"),
  value = c(paste(cfg$ranks, collapse = " "),
            sprintf("%d (%s)", cfg$collapse_level, cfg$da_rank),
            sprintf("%d (%s)", cfg$overview_level, cfg$overview_rank),
            if (is.null(cfg$rarefy_depth)) "smallest library" else cfg$rarefy_depth,
            cfg$seed, cfg$group,
            if (is.null(cfg$group_ref)) "first level alphabetically" else cfg$group_ref),
  source = unname(src))
write.csv(used, out("settings_used.csv"), row.names = FALSE)
message("Settings used by downstream_analysis.R:")
message(paste(sprintf("  %-26s = %-45s (%s)", used$setting, used$value, used$source), collapse = "\n"))
set.seed(cfg$seed)

# ---- 1. Import ---------------------------------------------------------------

read_qiime_metadata <- function(path) {
  md <- read.delim(path, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE,
                   colClasses = "character")
  md <- md[!grepl("^#q2:", md[[1]]), , drop = FALSE]   # drop the optional #q2:types row
  rownames(md) <- md[[1]]
  md <- md[, -1, drop = FALSE]
  # Convert columns that are entirely numeric (ignoring blanks) back to numbers.
  for (col in names(md)) {
    v <- md[[col]]
    num <- suppressWarnings(as.numeric(v))
    if (all(is.na(v) | v == "" | !is.na(num))) md[[col]] <- num
  }
  md
}

read_feature_table <- function(path) {
  # biom convert --to-tsv: line 1 is a comment, line 2 the header ("#OTU ID", samples...)
  tab <- read.delim(path, skip = 1, check.names = FALSE, comment.char = "",
                    row.names = 1, stringsAsFactors = FALSE)
  as.matrix(tab)
}

read_taxonomy <- function(path, ranks) {
  tx <- read.delim(path, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  split <- strsplit(tx$Taxon, ";\\s*")
  depth <- max(lengths(split))
  # Never truncate: a deeper taxonomy than `ranks` would shift every rank label [review #2].
  if (depth > length(ranks))
    stop(sprintf(paste("taxonomy.tsv has up to %d ranks but only %d rank names are configured",
                       "(%s). Set TAXONOMY_RANKS in pipeline/config.sh to match your reference."),
                 depth, length(ranks), paste(ranks, collapse = " ")))
  m <- t(vapply(split, function(x) {
    x <- trimws(x)
    # empty rank labels such as "g__", "sgr__", "dvs__" (any lowercase prefix)
    x[grepl("^[a-z]+__$", x) | x %in% c("", "Unassigned")] <- NA
    length(x) <- length(ranks)
    x
  }, character(length(ranks))))
  dimnames(m) <- list(tx[[1]], ranks)
  m   # a matrix keeps the rank names in tax_table() [E10]
}

counts   <- read_feature_table(file.path(cfg$export_dir, "feature-table.tsv"))
taxonomy <- read_taxonomy(file.path(cfg$export_dir, "taxonomy.tsv"), cfg$ranks)
metadata <- read_qiime_metadata(file.path(cfg$export_dir, "metadata.tsv"))

stopifnot(cfg$group %in% names(metadata))
shared <- intersect(colnames(counts), rownames(metadata))
if (length(setdiff(colnames(counts), shared)))
  say("Samples in table but not in metadata (dropped): %s",
      paste(setdiff(colnames(counts), shared), collapse = ", "))
metadata <- metadata[shared, , drop = FALSE]
metadata[[cfg$group]] <- factor(metadata[[cfg$group]])
if (!is.null(cfg$group_ref))
  metadata[[cfg$group]] <- relevel(metadata[[cfg$group]], ref = cfg$group_ref)

ps <- phyloseq(otu_table(counts[, shared, drop = FALSE], taxa_are_rows = TRUE),
               tax_table(taxonomy[rownames(counts), , drop = FALSE]),
               sample_data(metadata))

tree_path <- file.path(cfg$export_dir, "tree.nwk")
has_tree <- file.exists(tree_path)
if (has_tree) {
  tree <- ape::read.tree(tree_path)
  # SEPP trees also contain the reference sequences; keep only our ASVs.
  tree <- ape::keep.tip(tree, intersect(tree$tip.label, taxa_names(ps)))
  if (!ape::is.rooted(tree)) { say("Tree is unrooted; skipping UniFrac."); has_tree <- FALSE }
  else ps <- merge_phyloseq(ps, phy_tree(tree))
}
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
say("Imported %d samples x %d features (tree: %s)", nsamples(ps), ntaxa(ps), has_tree)

meta  <- data.frame(sample_data(ps), check.names = FALSE)
group <- meta[[cfg$group]]

# ---- 2. Depth QC and rarefaction curves -------------------------------------

depths <- sample_sums(ps)
write.csv(data.frame(sample = names(depths), reads = depths), out("sample_depths.csv"),
          row.names = FALSE)
depth <- if (is.null(cfg$rarefy_depth)) min(depths) else cfg$rarefy_depth
say("Library sizes: min %d, median %d, max %d; rarefaction depth %d",
    min(depths), round(median(depths)), max(depths), depth)

pdf(out("rarefaction_curves.pdf"), width = 9, height = 6)
rarecurve(t(as(otu_table(ps), "matrix")), step = max(1, floor(median(depths) / 50)),
          label = nsamples(ps) <= 40, cex = 0.5, xlab = "Reads", ylab = "ASVs")
abline(v = depth, lty = 2)
dev.off()

ps_even <- prune_samples(sample_sums(ps) >= depth, ps)

# ---- 3. Alpha diversity (averaged over repeated rarefaction) ----------------

alpha_once <- function(i) {
  r <- rarefy_even_depth(ps_even, sample.size = depth, rngseed = cfg$seed + i,
                         replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
  a <- estimate_richness(r, measures = c("Observed", "Shannon", "InvSimpson"))
  rownames(a) <- sample_names(r)  # estimate_richness mangles IDs via make.names()
  a$sample <- rownames(a)
  a
}
alpha <- bind_rows(lapply(seq_len(cfg$rarefy_iter), alpha_once)) |>
  group_by(sample) |>
  summarise(across(everything(), median), .groups = "drop") |>  # median, as q2-boots
  left_join(tibble::rownames_to_column(meta, "sample"), by = "sample")
write.csv(alpha, out("alpha_diversity.csv"), row.names = FALSE)

alpha_long <- pivot_longer(alpha, c(Observed, Shannon, InvSimpson),
                           names_to = "metric", values_to = "value")
group_test <- function(v, g) {
  g <- droplevels(factor(g))
  if (nlevels(g) == 2) wilcox.test(v ~ g, exact = FALSE)$p.value else kruskal.test(v, g)$p.value
}
alpha_tests <- alpha_long |>
  group_by(metric) |>
  summarise(test = if (nlevels(droplevels(.data[[cfg$group]])) == 2) "Wilcoxon" else "Kruskal-Wallis",
            p = group_test(value, .data[[cfg$group]]),
            .groups = "drop") |>
  mutate(q = p.adjust(p, "BH"))
write.csv(alpha_tests, out("alpha_tests.csv"), row.names = FALSE)

ggsave(out("alpha_diversity.pdf"), width = 9, height = 4,
       ggplot(alpha_long, aes(.data[[cfg$group]], value, fill = .data[[cfg$group]])) +
         geom_boxplot(outlier.shape = NA, alpha = 0.6) +
         geom_jitter(width = 0.15, size = 1) +
         facet_wrap(~ metric, scales = "free_y") +
         labs(x = NULL, y = NULL,
              subtitle = sprintf("Median of %d rarefactions to %d reads", cfg$rarefy_iter, depth)) +
         theme_bw() + theme(legend.position = "none"))

# ---- 4. Beta diversity, PERMANOVA and dispersion ----------------------------

rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_rare1 <- rarefy_even_depth(ps_even, sample.size = depth, rngseed = cfg$seed,
                              replace = FALSE, verbose = FALSE)
dists <- list(
  bray_curtis      = phyloseq::distance(rel, method = "bray"),
  robust_aitchison = vegdist(t(as(otu_table(ps), "matrix")), method = "robust.aitchison"),
  # presence/absence metrics depend on depth, so they use a rarefied table
  jaccard          = phyloseq::distance(ps_rare1, method = "jaccard", binary = TRUE)
)
if (has_tree) {
  dists$weighted_unifrac   <- UniFrac(rel, weighted = TRUE)
  dists$unweighted_unifrac <- UniFrac(ps_rare1, weighted = FALSE)
}

rhs <- paste(c(cfg$group, cfg$covariates), collapse = " + ")
beta_tests <- list()
for (nm in names(dists)) {
  d <- dists[[nm]]
  md <- meta[labels(d), , drop = FALSE]
  # by = "margin": each term is tested given all others, so term order does not matter
  # (QIIME 2 `diversity adonis` uses sequential terms instead) [E14].
  perm <- adonis2(as.formula(paste("d ~", rhs)), data = md, by = "margin",
                  permutations = cfg$permutations)
  disp <- permutest(betadisper(d, md[[cfg$group]]), permutations = cfg$permutations)
  beta_tests[[nm]] <- data.frame(
    distance = nm, term = rownames(perm)[seq_len(nrow(perm) - 2)],
    R2 = perm$R2[seq_len(nrow(perm) - 2)], p = perm$`Pr(>F)`[seq_len(nrow(perm) - 2)],
    group_dispersion_p = disp$tab$`Pr(>F)`[1])  # PERMDISP on `group`, same for every term

  pc <- cmdscale(d, k = 2, eig = TRUE)
  var_expl <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  pts <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], md, check.names = FALSE)
  p <- ggplot(pts, aes(PC1, PC2, colour = .data[[cfg$group]])) +
    geom_point(size = 2.5) +
    labs(x = sprintf("PCo1 (%.1f%%)", var_expl[1]), y = sprintf("PCo2 (%.1f%%)", var_expl[2]),
         title = gsub("_", " ", nm)) +
    theme_bw()
  # Ellipses only make sense for categorical groups with enough points [E13].
  if (all(table(pts[[cfg$group]]) >= 4)) p <- p + stat_ellipse(linetype = 2)
  ggsave(out(sprintf("pcoa_%s.pdf", nm)), p, width = 6, height = 5)
}
beta_tests <- bind_rows(beta_tests)
write.csv(beta_tests, out("beta_permanova_dispersion.csv"), row.names = FALSE)

# ---- 5. Composition barplots (top taxa + "Other") ---------------------------

plot_composition <- function(rank) {
  g <- tax_glom(rel, taxrank = rank, NArm = FALSE)
  df <- psmelt(g)
  df$taxon <- ifelse(is.na(df[[rank]]), "Unassigned", df[[rank]])
  top <- df |> group_by(taxon) |> summarise(m = mean(Abundance), .groups = "drop") |>
    filter(taxon != "Unassigned") |> slice_max(m, n = cfg$top_n) |> pull(taxon)
  df <- df |>
    mutate(taxon = ifelse(taxon %in% top | taxon == "Unassigned", taxon, "Other")) |>
    group_by(Sample, taxon, .data[[cfg$group]]) |>
    summarise(Abundance = sum(Abundance), .groups = "drop") |>
    mutate(taxon = factor(taxon, levels = c(top, "Other", "Unassigned")))
  # Pooled categories in neutral greys so they never share a hue with a real taxon.
  fills <- c(setNames(scales::hue_pal()(length(top)), top), Other = "grey75", Unassigned = "grey45")
  ggsave(out(sprintf("barplot_%s.pdf", tolower(rank))), width = 11, height = 5,
         ggplot(df, aes(Sample, Abundance, fill = taxon)) +
           geom_col(width = 0.95) +
           scale_fill_manual(values = fills) +
           facet_grid(~ .data[[cfg$group]], scales = "free_x", space = "free_x") +
           scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
           labs(x = NULL, y = "Relative abundance", fill = rank) +
           theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)))
}
invisible(lapply(unique(c(cfg$overview_rank, cfg$da_rank)), plot_composition))

# ---- 6. Differential abundance: ANCOM-BC2 + ALDEx2 + MaAsLin3 consensus -----

ps_da <- tax_glom(ps, taxrank = cfg$da_rank, NArm = FALSE)
lineage <- apply(as(tax_table(ps_da), "matrix")[, seq_len(match(cfg$da_rank, cfg$ranks)),
                                                drop = FALSE], 1,
                 function(x) paste(ifelse(is.na(x), "NA", x), collapse = ";"))
taxa_names(ps_da) <- make.unique(unname(lineage))
prev <- rowMeans(as(otu_table(ps_da), "matrix") > 0)
ps_da <- prune_taxa(prev >= cfg$da_prevalence, ps_da)
counts_da <- as(otu_table(ps_da), "matrix")            # unrarefied counts for every tool [E12]
meta_da <- data.frame(sample_data(ps_da), check.names = FALSE)
lv <- levels(meta_da[[cfg$group]])
two_groups <- length(lv) == 2
say("DA at %s level: %d taxa x %d samples; groups: %s",
    cfg$da_rank, ntaxa(ps_da), nsamples(ps_da), paste(lv, collapse = " vs "))

da <- list()

# ALDEx2 (CLR; Wilcoxon for two groups, Kruskal-Wallis otherwise) [MH_2026.3-R]
if (requireNamespace("ALDEx2", quietly = TRUE)) {
  clr <- ALDEx2::aldex.clr(counts_da, as.character(meta_da[[cfg$group]]),
                           mc.samples = 128, denom = "all", verbose = FALSE)
  if (two_groups) {
    tt <- ALDEx2::aldex.ttest(clr)
    ef <- ALDEx2::aldex.effect(clr, verbose = FALSE)
    # ALDEx2 orders the condition labels alphabetically: effect = (2nd - 1st) of sort(lv).
    # ANCOM-BC2 and MaAsLin3 report (other level - reference), with lv[1] the reference,
    # so flip ALDEx2's sign when the reference is not alphabetically first.
    aldex_sign <- if (identical(lv[1], sort(lv)[1])) 1 else -1
    da$ALDEx2 <- data.frame(feature = rownames(tt), q = tt$wi.eBH,
                            effect = aldex_sign * ef$effect)
  } else {
    kw <- ALDEx2::aldex.kw(clr)
    da$ALDEx2 <- data.frame(feature = rownames(kw), q = kw$kw.eBH, effect = NA_real_)
  }
}

# MaAsLin3 (TSS + log; read depth as covariate, as in [MH_2026.3-R]) [E8]
if (requireNamespace("maaslin3", quietly = TRUE)) {
  mm <- meta_da
  mm$log_depth <- log(colSums(counts_da)[rownames(mm)])
  fit <- maaslin3::maaslin3(
    input_data = as.data.frame(t(counts_da)), input_metadata = mm,
    output = out("maaslin3"),
    formula = paste("~", paste(c(cfg$group, cfg$covariates, "log_depth"), collapse = " + ")),
    normalization = "TSS", transform = "LOG", augment = TRUE, standardize = TRUE,
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
    plot_summary_plot = FALSE, plot_associations = FALSE, verbosity = "ERROR", cores = 1)
  res <- fit$fit_data_abundance$results
  res <- res[res$metadata == cfg$group & is.na(res$error), ]
  # one row per feature: the most significant non-reference level
  res <- res[order(res$qval_individual), ]
  res <- res[!duplicated(res$feature), ]
  da$MaAsLin3 <- data.frame(feature = res$feature, q = res$qval_individual, effect = res$coef)
}

# ANCOM-BC2 (bias-corrected log-linear model) [gut_to_soil_tutorial]
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  # Matrix input (with meta_data) avoids the extra `microbiome` dependency of phyloseq input.
  ab <- ANCOMBC::ancombc2(
    data = counts_da, taxa_are_rows = TRUE, meta_data = meta_da,
    fix_formula = paste(c(cfg$group, cfg$covariates), collapse = " + "),
    p_adj_method = "BH", prv_cut = 0, lib_cut = 0, group = cfg$group,
    struc_zero = FALSE, neg_lb = FALSE, alpha = cfg$da_q, pseudo_sens = TRUE,
    global = FALSE, n_cl = 1, verbose = FALSE)$res
  qcols <- grep(paste0("^q_", cfg$group), names(ab), value = TRUE)
  lcols <- sub("^q_", "lfc_", qcols)
  best <- apply(as.matrix(ab[, qcols, drop = FALSE]), 1, which.min)
  qbest <- as.matrix(ab[, qcols, drop = FALSE])[cbind(seq_len(nrow(ab)), best)]
  # Require the pseudo-count sensitivity analysis to pass, when it was run.
  sscols <- sub("^q_", "passed_ss_", qcols)
  if (all(sscols %in% names(ab))) {
    ss <- as.matrix(ab[, sscols, drop = FALSE])[cbind(seq_len(nrow(ab)), best)]
    qbest[!ss] <- 1
  }
  da$ANCOMBC2 <- data.frame(
    feature = ab$taxon, q = qbest,
    effect = as.matrix(ab[, lcols, drop = FALSE])[cbind(seq_len(nrow(ab)), best)])
}

if (length(da) == 0) stop("No DA package installed (ALDEx2, maaslin3, ANCOMBC).")
da_all <- bind_rows(lapply(names(da), function(m) mutate(da[[m]], method = m)))
write.csv(da_all, out("da_all_methods.csv"), row.names = FALSE)

consensus <- da_all |>
  mutate(sig = !is.na(q) & q <= cfg$da_q, direction = sign(effect)) |>
  group_by(feature) |>
  summarise(n_methods_sig = sum(sig),
            methods = paste(method[sig], collapse = ","),
            # agreement on direction among significant calls (two-group designs)
            direction_agrees = if (two_groups) n_distinct(direction[sig & !is.na(direction)]) <= 1
                               else NA,
            .groups = "drop") |>
  filter(n_methods_sig >= min(cfg$da_min_methods, length(da))) |>
  filter(is.na(direction_agrees) | direction_agrees) |>
  arrange(desc(n_methods_sig))
write.csv(consensus, out("da_consensus.csv"), row.names = FALSE)
say("DA consensus (>= %d of %d methods, q <= %.2f): %d taxa",
    min(cfg$da_min_methods, length(da)), length(da), cfg$da_q, nrow(consensus))

if (nrow(consensus) > 0) {
  rel_da <- transform_sample_counts(ps_da, function(x) x / sum(x))
  show <- head(consensus$feature, 12)
  df <- psmelt(prune_taxa(show, rel_da))
  df$label <- sub(".*;", "", df$OTU)
  ggsave(out("da_consensus_boxplots.pdf"), width = 10, height = 7,
         ggplot(df, aes(.data[[cfg$group]], Abundance, fill = .data[[cfg$group]])) +
           geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 1) +
           facet_wrap(~ label, scales = "free_y") +
           scale_y_continuous(labels = scales::percent) +
           labs(x = NULL, y = "Relative abundance") +
           theme_bw() + theme(legend.position = "none"))
}

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
say("Done. Results in %s", cfg$out_dir)
