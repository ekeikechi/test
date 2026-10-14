#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Design-aware statistics for nested / split-plot field experiments (stage 11b).
#
# R/downstream_analysis.R compares groups of independent samples. Field designs where
# treatments are applied to plots, sub-treatments to subplots within plots and samples
# are taken at several depths in each subplot need the plots and subplots in the model,
# otherwise samples from one plot count as independent replicates (pseudo-replication).
# This script reads one design file per analysis (R/designs/*.tsv) and runs:
#
#   1. Sample subset (e.g. depth_cm == 0), design checks, the stratum of every term
#   2. Depth rules: MIN_SAMPLE_DEPTH, main depth = smallest remaining sample, and a
#      sensitivity depth keeping >= DEPTH_SENSITIVITY_MIN_KEEP of every group
#      (same rule as pipeline/choose_depth.py)
#   3. Alpha diversity: median of repeated rarefactions; linear mixed model with random
#      plot and subplot intercepts (lme4/lmerTest), type III tests, Satterthwaite df
#   4. Beta diversity: hierarchical PERMANOVA. Each term is tested in its stratum, on
#      units of that stratum, with permutations restricted accordingly:
#        plot stratum     plot-mean profiles,    d ~ plot terms,              free
#        subplot stratum  subplot-mean profiles, d ~ plot + subplot terms,    within plots
#        sample stratum   samples,               d ~ subplot + sample terms,  within subplots
#      plus PERMDISP for every factor in its stratum
#   5. Differential abundance per model coefficient: ANCOM-BC2 (rand_formula) and
#      MaAsLin3 (random effects) are primary; consensus = q <= DA_Q in both with the same
#      sign. ALDEx2 (no random effects possible) is reported as supporting evidence only.
#   6. Alpha diversity and Jaccard PERMANOVA again at the sensitivity depth
#
# Usage : Rscript R/design_analysis.R DESIGN_FILE [OUT_DIR]
# Inputs: a stage 10-style export (feature-table.tsv, taxonomy.tsv, metadata.tsv,
#         analysis_settings.tsv), e.g. results/10_export or an R/pr2_traits.R subset.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
})
say <- function(...) message(sprintf(...))
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript R/design_analysis.R DESIGN_FILE [OUT_DIR]")

# ---- 0. Design file and shared settings --------------------------------------

read_kv <- function(path) {
  x <- read.delim(path, colClasses = "character", quote = "", comment.char = "#",
                  na.strings = character(0), strip.white = TRUE)
  stopifnot(identical(names(x)[1:2], c("setting", "value")))
  setNames(x$value, x$setting)
}
design <- read_kv(args[[1]])
known <- c("name", "description", "export_dir", "out_dir", "sample_filter", "fixed", "da_fixed",
           "plot", "subplot", "reference", "expect_strata", "sampling_depth", "min_sample_depth",
           "depth_groups", "depth_min_keep", "da_rank", "da_prevalence", "da_q", "permutations",
           "rarefy_iter", "distances", "da_methods", "cores", "seed")
if (length(bad <- setdiff(names(design), known)))
  stop(sprintf("%s: unknown setting(s): %s", args[[1]], paste(bad, collapse = ", ")))

value <- list(); source_of <- list()
setv <- function(key, val, src) { value[[key]] <<- val; source_of[[key]] <<- src }
dget <- function(k) if (k %in% names(design) && nzchar(design[[k]])) design[[k]] else NULL
for (k in names(design)) if (nzchar(design[[k]])) setv(k, design[[k]], "design file")
if (is.null(dget("name"))) setv("name", sub("\\.tsv$", "", basename(args[[1]])), "file name")
if (is.null(dget("export_dir"))) setv("export_dir", "results/10_export", "default")
if (is.null(dget("fixed"))) stop("the design file must give `fixed` (the fixed-effect formula)")
if (is.null(dget("plot"))) stop("the design file must give `plot` (the whole-plot unit column)")

settings_path <- file.path(value$export_dir, "analysis_settings.tsv")
st <- if (file.exists(settings_path)) read_kv(settings_path) else character(0)
sget <- function(k) if (k %in% names(st) && nzchar(st[[k]])) st[[k]] else NULL
inherit <- function(key, setting, default) {  # design file > analysis_settings.tsv > default
  if (!is.null(dget(key))) return(invisible())
  if (!is.null(v <- sget(setting))) setv(key, v, "analysis_settings.tsv")
  else setv(key, default, "default")
}
inherit("sampling_depth", "sampling_depth", "auto")
inherit("min_sample_depth", "min_sample_depth", "")
inherit("depth_groups", "depth_sensitivity_groups", "")
inherit("depth_min_keep", "depth_sensitivity_min_keep", "0.8")
inherit("reference", "da_reference", "")
inherit("seed", "random_seed", "711")
ranks <- if (!is.null(v <- sget("taxonomy_ranks"))) strsplit(v, "\\s+")[[1]] else
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (is.null(dget("da_rank")))
  setv("da_rank", ranks[if (!is.null(v <- sget("collapse_level"))) as.integer(v) else 6],
       if (!is.null(sget("collapse_level"))) "analysis_settings.tsv" else "default")
defaults <- c(da_fixed = value$fixed, da_prevalence = "0.1", da_q = "0.1", permutations = "999",
              rarefy_iter = "100", distances = "bray_curtis robust_aitchison jaccard",
              da_methods = "ANCOMBC2 MaAsLin3 ALDEx2", cores = "1", sample_filter = "",
              subplot = "", expect_strata = "", out_dir = "")
for (k in names(defaults)) if (is.null(dget(k))) setv(k, defaults[[k]], "default")
if (!nzchar(value$out_dir)) setv("out_dir", file.path("results/11_design", value$name), "default")
if (length(args) >= 2) setv("out_dir", args[[2]], "command line")

num <- function(k) as.numeric(value[[k]])
words <- function(k) strsplit(trimws(value[[k]]), "\\s+")[[1]]
dir.create(value$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(...) file.path(value$out_dir, ...)
write_tsv <- function(x, f) write.table(x, out(f), sep = "\t", quote = FALSE, row.names = FALSE,
                                        na = "")
seed <- as.integer(value$seed)
set.seed(seed)

# ---- 1. Import, subset, design checks ---------------------------------------

read_metadata <- function(path) {
  md <- read.delim(path, check.names = FALSE, comment.char = "", colClasses = "character")
  md <- md[!grepl("^#q2:", md[[1]]), , drop = FALSE]
  rownames(md) <- md[[1]]
  md <- md[, -1, drop = FALSE]
  for (col in names(md)) {
    v <- md[[col]]; n <- suppressWarnings(as.numeric(v))
    if (all(is.na(v) | v == "" | !is.na(n))) md[[col]] <- n
  }
  md
}
counts <- as.matrix(read.delim(file.path(value$export_dir, "feature-table.tsv"), skip = 1,
                               check.names = FALSE, comment.char = "", row.names = 1))
tax <- read.delim(file.path(value$export_dir, "taxonomy.tsv"), check.names = FALSE,
                  colClasses = "character")
meta <- read_metadata(file.path(value$export_dir, "metadata.tsv"))

fixed_f <- as.formula(paste("~", value$fixed))
da_f <- as.formula(paste("~", value$da_fixed))
vars <- unique(c(all.vars(fixed_f), all.vars(da_f)))
plot_col <- value$plot
sub_col <- if (nzchar(value$subplot)) value$subplot else NULL
need <- c(vars, plot_col, sub_col)
if (length(miss <- setdiff(need, names(meta))))
  stop(sprintf("metadata.tsv has no column(s): %s", paste(miss, collapse = ", ")))

shared <- intersect(colnames(counts), rownames(meta))
if (length(d <- setdiff(colnames(counts), shared)))
  say("Samples in the table but not in the metadata (dropped): %s", paste(d, collapse = ", "))
meta <- meta[shared, , drop = FALSE]
n_all <- nrow(meta)
if (nzchar(value$sample_filter)) {
  keep <- eval(parse(text = value$sample_filter), meta, baseenv())
  if (!is.logical(keep) || length(keep) != nrow(meta))
    stop(sprintf("sample_filter `%s` must give TRUE/FALSE per sample", value$sample_filter))
  meta <- meta[!is.na(keep) & keep, , drop = FALSE]
}
say("%s: %d of %d samples pass the sample filter `%s`", value$name, nrow(meta), n_all,
    if (nzchar(value$sample_filter)) value$sample_filter else "(none)")
if (nrow(meta) < 4) stop("fewer than 4 samples left after the sample filter")
counts <- counts[, rownames(meta), drop = FALSE]

# Factors: character columns used in a formula; levels in natural order ("1cm" < "10cm"),
# then the references from DA_REFERENCE / the design file ("column::level").
natural_levels <- function(x) {
  u <- unique(as.character(x))
  lead <- suppressWarnings(as.numeric(sub("^(-?[0-9.]+).*$", "\\1", u)))
  if (all(!is.na(lead))) u[order(lead, u)] else sort(u)
}
refs <- if (nzchar(value$reference)) strsplit(value$reference, "\\s+")[[1]] else character(0)
for (v in vars) if (!is.numeric(meta[[v]])) {
  meta[[v]] <- factor(meta[[v]], levels = natural_levels(meta[[v]]))
  hit <- grep(paste0("^", v, "::"), refs, value = TRUE)
  if (length(hit)) {
    lv <- sub("^[^:]*::", "", hit[1])
    if (!lv %in% levels(meta[[v]]))
      stop(sprintf("reference level `%s` is not a level of %s (%s)", lv, v,
                   paste(levels(meta[[v]]), collapse = ", ")))
    meta[[v]] <- relevel(meta[[v]], ref = lv)
  }
  meta[[v]] <- droplevels(meta[[v]])
}
meta[[plot_col]] <- factor(meta[[plot_col]])
if (!is.null(sub_col)) {
  meta[[sub_col]] <- factor(meta[[sub_col]])
  per_sub <- tapply(as.character(meta[[plot_col]]), meta[[sub_col]], function(x) length(unique(x)))
  if (any(per_sub > 1))
    stop(sprintf(paste("subplot IDs must be unique across plots; %s occurs in several plots.",
                       "Use e.g. plot_exclosure as the subplot ID."),
                 paste(names(per_sub)[per_sub > 1], collapse = ", ")))
}

# Stratum of each variable: constant within plots -> plot; within subplots -> subplot;
# otherwise sample. A term's stratum is the lowest stratum among its variables.
levels_order <- c("plot", "subplot", "sample")
unit_col <- c(plot = plot_col, subplot = if (is.null(sub_col)) NA else sub_col, sample = NA)
one_per_subplot <- !is.null(sub_col) && all(table(meta[[sub_col]]) == 1)
constant_within <- function(v, g) all(tapply(as.character(meta[[v]]), meta[[g]],
                                             function(x) length(unique(x))) <= 1)
var_stratum <- vapply(vars, function(v) {
  if (constant_within(v, plot_col)) "plot"
  else if (!is.null(sub_col) && constant_within(v, sub_col)) "subplot"
  else "sample"
}, character(1))
# with one sample per subplot (e.g. a single depth) the subplot and sample strata coincide
if (one_per_subplot) var_stratum[var_stratum == "sample"] <- "subplot"
term_stratum <- function(f) {
  labs <- attr(terms(f), "term.labels")
  setNames(vapply(labs, function(t) {
    vs <- all.vars(as.formula(paste("~", t)))
    levels_order[max(match(var_stratum[vs], levels_order))]
  }, character(1)), labs)
}
strata <- term_stratum(fixed_f)
if (nzchar(value$expect_strata)) {
  for (item in strsplit(value$expect_strata, "\\s+")[[1]]) {
    kv <- strsplit(item, "=", fixed = TRUE)[[1]]
    if (!kv[1] %in% names(var_stratum)) next   # e.g. depth when only 0 cm is analysed
    got <- var_stratum[[kv[1]]]
    if (got != kv[2] && !(one_per_subplot && kv[2] == "sample" && got == "subplot"))
      stop(sprintf(paste("%s varies at the %s level, but the design file expects %s.",
                         "Check the plot/subplot IDs in metadata.tsv."), kv[1], got, kv[2]))
  }
}
unit_label <- function(s) switch(s, plot = sprintf("%s means", plot_col),
                                 subplot = if (one_per_subplot) "samples (one per subplot)"
                                           else sprintf("%s means", sub_col),
                                 sample = "samples")
write_tsv(data.frame(term = names(strata), stratum = unname(strata),
                     tested_on = vapply(strata, unit_label, ""),
                     permutations = vapply(strata, function(s) switch(s,
                       plot = "free", subplot = sprintf("within %s", plot_col),
                       sample = sprintf("within %s", sub_col %||% plot_col)), "")),
          "design_strata.tsv")
cells <- as.data.frame(table(meta[, all.vars(fixed_f), drop = FALSE]))
write_tsv(cells, "design_cells.tsv")
if (length(unique(cells$Freq)) > 1)
  say("Note: the design is unbalanced (%d-%d samples per cell); see design_cells.tsv.",
      min(cells$Freq), max(cells$Freq))
say("Strata: %s", paste(sprintf("%s=%s", names(strata), strata), collapse = ", "))

# ---- 2. Depth rules ------------------------------------------------------------

depths <- colSums(counts)
removed <- character(0)
if (grepl("^[0-9]+$", value$min_sample_depth) && any(depths < as.numeric(value$min_sample_depth))) {
  removed <- names(depths)[depths < as.numeric(value$min_sample_depth)]
  say("Removing %d sample(s) below MIN_SAMPLE_DEPTH = %s reads: %s", length(removed),
      value$min_sample_depth, paste(removed, collapse = ", "))
}
write_tsv(data.frame(sample = removed, reads = depths[removed]), "samples_removed.tsv")
meta <- meta[setdiff(rownames(meta), removed), , drop = FALSE]
for (v in c(vars, plot_col, sub_col)) if (is.factor(meta[[v]])) meta[[v]] <- droplevels(meta[[v]])
counts <- counts[, rownames(meta), drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
depths <- colSums(counts)

# Same rule as pipeline/choose_depth.py
choose_depths <- function(depths, meta, groups, min_keep) {
  main <- as.integer(min(depths))
  if (!length(groups)) return(list(main = main, sens = NA_integer_, by_group = NULL))
  g <- interaction(meta[names(depths), groups, drop = FALSE], sep = " / ", drop = TRUE)
  allowed <- tapply(depths, g, function(ds) {
    ds <- sort(ds, decreasing = TRUE); ds[max(1, ceiling(min_keep * length(ds)))]
  })
  sens <- as.integer(min(allowed))
  by_group <- data.frame(group = names(allowed), samples = as.vector(table(g)[names(allowed)]),
                         kept_at_main = as.vector(tapply(depths >= main, g, sum)[names(allowed)]),
                         kept_at_sensitivity = as.vector(tapply(depths >= sens, g, sum)[names(allowed)]),
                         dropped_at_sensitivity = vapply(names(allowed), function(x)
                           paste(sort(names(depths)[g == x & depths < sens]), collapse = ","), ""))
  list(main = main, sens = sens, by_group = by_group,
       limiting = paste(sort(names(allowed)[as.integer(allowed) == sens]), collapse = "; "))
}
groups <- if (nzchar(value$depth_groups)) words("depth_groups") else character(0)
groups <- groups[groups %in% names(meta)]
dc <- choose_depths(depths, meta, groups, num("depth_min_keep"))
main_depth <- if (grepl("^[0-9]+$", value$sampling_depth)) as.integer(value$sampling_depth) else dc$main
sens_depth <- if (!is.na(dc$sens) && dc$sens != main_depth) dc$sens else NA_integer_
write_tsv(data.frame(setting = c("samples", "main_depth", "main_rule", "sensitivity_depth",
                                 "sensitivity_rule", "limiting_groups",
                                 "samples_kept_at_sensitivity"),
                     value = c(length(depths), main_depth,
                               if (grepl("^[0-9]+$", value$sampling_depth)) "SAMPLING_DEPTH" else
                                 "smallest remaining sample (all samples kept)",
                               if (is.na(dc$sens)) "" else dc$sens,
                               if (!length(groups)) "" else sprintf(
                                 "deepest depth keeping >= %s of samples in every group of: %s",
                                 value$depth_min_keep, paste(groups, collapse = " x ")),
                               dc$limiting %||% "",
                               if (is.na(dc$sens)) "" else sum(depths >= dc$sens))),
          "depth_choice.tsv")
if (!is.null(dc$by_group)) write_tsv(dc$by_group, "depth_by_group.tsv")
say("%d samples x %d features; depths: main %d, sensitivity %s", ncol(counts), nrow(counts),
    main_depth, if (is.na(dc$sens)) "none (no depth groups)" else dc$sens)

# ---- Random-effect structure ----------------------------------------------------

# A random intercept needs more than one sample per level and fewer levels than samples.
usable <- function(col) {
  if (is.null(col)) return(FALSE)
  n <- table(meta[[col]]); nlevels(droplevels(meta[[col]])) > 1 && any(n > 1) &&
    nlevels(droplevels(meta[[col]])) < nrow(meta)
}
random_cols <- Filter(usable, c(plot_col, sub_col))
random_terms <- if (length(random_cols)) paste(sprintf("(1 | %s)", random_cols), collapse = " + ") else ""
say("Random intercepts: %s", if (nzchar(random_terms)) random_terms else "none")

# ---- 3. Alpha diversity ------------------------------------------------------------

rarefied_alpha <- function(depth, iter) {
  keep <- colnames(counts)[depths >= depth]
  x <- t(counts[, keep, drop = FALSE])
  runs <- lapply(seq_len(iter), function(i) {
    set.seed(seed + i)
    r <- rrarefy(x, depth)
    cbind(Observed = specnumber(r), Shannon = diversity(r, "shannon"),
          InvSimpson = diversity(r, "invsimpson"))
  })
  med <- apply(simplify2array(runs), c(1, 2), median)   # median, as q2-boots
  data.frame(sample = keep, med, check.names = FALSE)
}
transforms <- c(Observed = "log", Shannon = "identity", InvSimpson = "log")

fit_alpha <- function(alpha, suffix) {
  d <- cbind(alpha, meta[alpha$sample, , drop = FALSE])
  op <- options(contrasts = c("contr.sum", "contr.poly"))   # type III tests need sum-to-zero
  on.exit(options(op))
  have_lmm <- requireNamespace("lmerTest", quietly = TRUE)
  rc <- Filter(function(col) {
    n <- table(droplevels(d[[col]])); length(n) > 1 && any(n > 1) && length(n) < nrow(d)
  }, random_cols)
  res <- list(); vc <- list()
  for (m in names(transforms)) {
    d$y <- if (transforms[[m]] == "log") log(d[[m]]) else d[[m]]
    rhs <- paste(c(deparse1(fixed_f[[2]]), sprintf("(1 | %s)", rc)), collapse = " + ")
    if (length(rc) && have_lmm) {
      fit <- lmerTest::lmer(as.formula(paste("y ~", rhs)), data = d, REML = TRUE)
      a <- anova(fit, type = 3, ddf = "Satterthwaite")
      res[[m]] <- data.frame(metric = m, transform = transforms[[m]], term = rownames(a),
                             NumDF = a$NumDF, DenDF = a$DenDF, F = a$`F value`,
                             p = a$`Pr(>F)`, singular = lme4::isSingular(fit))
      v <- as.data.frame(lme4::VarCorr(fit))
      vc[[m]] <- data.frame(metric = m, component = v$grp, variance = v$vcov,
                            share = v$vcov / sum(v$vcov))
    } else {
      if (length(rc)) say("lmerTest not installed: alpha diversity uses lm WITHOUT plot terms")
      fit <- lm(as.formula(paste("y ~", deparse1(fixed_f[[2]]))), data = d)
      a <- drop1(fit, scope = ~ ., test = "F")[-1, ]
      res[[m]] <- data.frame(metric = m, transform = transforms[[m]], term = rownames(a),
                             NumDF = a$Df, DenDF = fit$df.residual, F = a$`F value`,
                             p = a$`Pr(>F)`, singular = NA)
    }
  }
  res <- do.call(rbind, res)
  res$stratum <- unname(strata[res$term])
  write_tsv(res, sprintf("alpha_lmm%s.tsv", suffix))
  if (length(vc)) write_tsv(do.call(rbind, vc), sprintf("alpha_variance_components%s.tsv", suffix))
  res
}

iter <- as.integer(value$rarefy_iter)
alpha <- rarefied_alpha(main_depth, iter)
write_tsv(merge(alpha, cbind(sample = rownames(meta), meta), by = "sample"), "alpha_diversity.tsv")
alpha_res <- fit_alpha(alpha, "")

# Figure: metrics by the treatment combination of the plot/subplot factors, per sample factor
fac_vars <- names(var_stratum)[vapply(names(var_stratum), function(v) is.factor(meta[[v]]), TRUE)]
x_vars <- fac_vars[var_stratum[fac_vars] != "sample"]
facet_vars <- fac_vars[var_stratum[fac_vars] == "sample"]
along <- reshape(merge(alpha, cbind(sample = rownames(meta), meta), by = "sample"),
                 direction = "long", varying = names(transforms), v.names = "value",
                 timevar = "metric", times = names(transforms))
along$x <- if (length(x_vars)) interaction(along[, x_vars, drop = FALSE], sep = "\n") else "all"
p <- ggplot(along, aes(x, value)) + geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  labs(x = NULL, y = NULL, subtitle = sprintf("Median of %d rarefactions to %d reads", iter, main_depth)) +
  theme_bw() + theme(axis.text.x = element_text(size = 6))
p <- p + if (length(facet_vars)) facet_grid(as.formula(paste("metric ~", paste(facet_vars, collapse = "+"))),
                                            scales = "free_y") else facet_wrap(~ metric, scales = "free_y")
ggsave(out("alpha_diversity.pdf"), p, width = 11, height = 7)

# ---- 4. Hierarchical PERMANOVA ------------------------------------------------------

prop <- sweep(counts, 2, depths, "/")
rare_presence <- function(depth) {  # 0/1 after one rarefaction (taxa x samples)
  keep <- colnames(counts)[depths >= depth]
  set.seed(seed)
  (t(rrarefy(t(counts[, keep, drop = FALSE]), depth)) > 0) * 1
}
# Profiles are taxa x units. Bray-Curtis and robust Aitchison use (mean) proportions.
# Jaccard uses presence/absence: binary Jaccard for samples; for plot or subplot means the
# profile is the share of samples where each taxon occurs, compared with quantitative
# Jaccard (2B / (1 + B)), which equals binary Jaccard on 0/1 data. (Binary Jaccard on
# pooled samples would only compare which taxa occur anywhere in a plot.)
dist_of <- function(name, profile) {
  x <- t(profile)
  switch(name,
         bray_curtis = vegdist(x, "bray"),
         robust_aitchison = vegdist(x, "robust.aitchison"),
         jaccard = vegdist(x, "jaccard"),
         stop("unknown distance ", name))
}
term_key <- function(t) vapply(strsplit(t, ":", fixed = TRUE),
                               function(v) paste(sort(v), collapse = ":"), "")
aggregate_units <- function(profile, unit) {  # mean profile per unit (taxa x units)
  if (is.null(unit)) return(profile)
  u <- as.character(meta[colnames(profile), unit])
  agg <- rowsum(t(profile), u)
  t(agg / as.vector(table(u)[rownames(agg)]))
}
unit_meta <- function(units, unit) {
  if (is.null(unit)) return(meta[units, , drop = FALSE])
  m <- meta[match(units, as.character(meta[[unit]])), , drop = FALSE]
  rownames(m) <- units
  m
}
level_spec <- function(s) {  # aggregation unit, conditioning factor, permutation blocks
  switch(s,
         plot = list(unit = plot_col, cond = NULL),
         subplot = list(unit = if (one_per_subplot) NULL else sub_col, cond = plot_col),
         sample = list(unit = NULL, cond = sub_col %||% plot_col))
}

permanova <- function(dist_names, depth, suffix) {
  out_rows <- list(); disp_rows <- list()
  nperm <- as.integer(value$permutations)
  for (dn in dist_names) {
    base <- if (dn == "jaccard") rare_presence(depth) else prop
    for (s in levels_order) {
      tl <- names(strata)[strata == s]
      if (!length(tl)) next
      sp <- level_spec(s)
      prof <- aggregate_units(base, sp$unit)
      md <- unit_meta(colnames(prof), sp$unit)
      d <- dist_of(dn, prof)
      rhs <- paste(c(sp$cond, tl), collapse = " + ")
      perms <- if (is.null(sp$cond)) how(nperm = nperm) else how(blocks = md[[sp$cond]], nperm = nperm)
      set.seed(seed)
      # sequential terms: the conditioning factor first, then the terms in formula order
      fit <- adonis2(as.formula(paste("d ~", rhs)), data = md, by = "terms",
                     permutations = perms)
      # R may reorder the variables of an interaction ("exclosure:nutrient")
      rows <- rownames(fit)
      keep <- term_key(rows) %in% term_key(tl)
      ss_resid <- fit$SumOfSqs[rows == "Residual"]
      out_rows[[length(out_rows) + 1]] <- data.frame(
        distance = dn, stratum = s, units = unit_label(s), n_units = nrow(md),
        conditioned_on = sp$cond %||% "",
        term = tl[match(term_key(rows[keep]), term_key(tl))], Df = fit$Df[keep],
        SumOfSqs = fit$SumOfSqs[keep], R2_of_stratum = fit$R2[keep], F = fit$F[keep],
        p = fit$`Pr(>F)`[keep], residual_df = fit$Df[rows == "Residual"],
        residual_SS = ss_resid)
      # PERMDISP for each factor that has its main effect in this stratum
      for (v in intersect(tl, names(var_stratum)[vapply(names(var_stratum),
                                                         function(x) is.factor(meta[[x]]), TRUE)])) {
        set.seed(seed)
        bd <- permutest(betadisper(d, md[[v]]), permutations = perms)
        disp_rows[[length(disp_rows) + 1]] <- data.frame(distance = dn, stratum = s, factor = v,
                                                         F = bd$tab$F[1], p = bd$tab$`Pr(>F)`[1])
      }
    }
  }
  res <- do.call(rbind, out_rows)
  write_tsv(res, sprintf("permanova_hierarchical%s.tsv", suffix))
  if (length(disp_rows)) write_tsv(do.call(rbind, disp_rows), sprintf("permdisp%s.tsv", suffix))
  res
}
dist_names <- words("distances")
perm_res <- permanova(dist_names, main_depth, "")

# PCoA of the samples for each distance, one page per factor
for (dn in dist_names) {
  d <- dist_of(dn, if (dn == "jaccard") rare_presence(main_depth) else prop)
  pc <- suppressWarnings(cmdscale(d, k = 2, eig = TRUE))
  if (is.null(dim(pc$points)) || ncol(pc$points) < 2) { say("PCoA of %s skipped", dn); next }
  ve <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  pts <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], meta[labels(d), , drop = FALSE],
                    check.names = FALSE)
  pdf(out(sprintf("pcoa_%s.pdf", dn)), width = 6.5, height = 5)
  for (v in fac_vars)
    print(ggplot(pts, aes(PC1, PC2, colour = .data[[v]])) + geom_point(size = 1.8, alpha = 0.8) +
            labs(x = sprintf("PCo1 (%.1f%%)", ve[1]), y = sprintf("PCo2 (%.1f%%)", ve[2]),
                 title = sprintf("%s, coloured by %s", gsub("_", " ", dn), v)) + theme_bw())
  dev.off()
}

# ---- 5. Differential abundance per coefficient -----------------------------------------

rank_i <- match(value$da_rank, ranks)
if (is.na(rank_i)) stop(sprintf("da_rank %s is not one of: %s", value$da_rank, paste(ranks, collapse = " ")))
lin <- vapply(strsplit(tax$Taxon, ";\\s*"), function(x) {
  x <- trimws(x); x[grepl("^[a-z]+__$", x) | x %in% c("", "Unassigned")] <- NA
  length(x) <- rank_i
  # unassigned ranks keep the deepest name, so they do not pool across lineages
  paste(ifelse(is.na(x), "NA", x), collapse = ";")
}, "")
names(lin) <- tax[[1]]
taxon <- lin[rownames(counts)]
if (anyNA(taxon)) stop("features without taxonomy: ", paste(head(rownames(counts)[is.na(taxon)]), collapse = ", "))
counts_da <- rowsum(counts, taxon)
prev <- rowMeans(counts_da > 0)
counts_da <- counts_da[prev >= num("da_prevalence"), , drop = FALSE]
meta_da <- meta[colnames(counts_da), , drop = FALSE]
mm <- model.matrix(da_f, meta_da)
coefs <- setdiff(colnames(mm), "(Intercept)")
say("DA at %s level: %d taxa x %d samples; model %s%s", value$da_rank, nrow(counts_da),
    ncol(counts_da), deparse1(da_f), if (nzchar(random_terms)) paste(" +", random_terms) else "")

# What each coefficient means (treatment coding, references as set above)
da_vars <- all.vars(da_f)
da_labels <- attr(terms(da_f), "term.labels")
coef_var <- function(p) {  # the variable a model-matrix column part belongs to
  da_vars[vapply(da_vars, function(x) if (is.numeric(meta_da[[x]])) p == x else
    startsWith(p, x) && substring(p, nchar(x) + 1) %in% levels(meta_da[[x]]), TRUE)][1]
}
coef_meaning <- vapply(coefs, function(cf) {
  parts <- strsplit(cf, ":", fixed = TRUE)[[1]]
  pv <- vapply(parts, coef_var, "")
  desc <- vapply(seq_along(parts), function(i) {
    v <- pv[[i]]
    if (is.numeric(meta_da[[v]])) sprintf("per unit of %s", v)
    else sprintf("%s: %s vs %s", v, substring(parts[[i]], nchar(v) + 1), levels(meta_da[[v]])[1])
  }, "")
  if (length(parts) == 1) {
    # main effects are simple effects at the reference level of interacting factors
    inter <- da_labels[grepl(":", da_labels) & vapply(da_labels, function(l)
      pv %in% strsplit(l, ":", fixed = TRUE)[[1]], TRUE)]
    others <- setdiff(unique(unlist(strsplit(inter, ":", fixed = TRUE))), pv)
    others <- others[!vapply(others, function(o) is.numeric(meta_da[[o]]), TRUE)]
    if (length(others)) desc <- sprintf("%s, at %s", desc, paste(sprintf("%s = %s", others,
      vapply(others, function(o) levels(meta_da[[o]])[1], "")), collapse = ", "))
    desc
  } else sprintf("interaction (difference of effects): %s", paste(desc, collapse = " x "))
}, "")
write_tsv(data.frame(coefficient = coefs, meaning = unname(coef_meaning)), "da_coefficients.tsv")

methods <- words("da_methods")
cores <- as.integer(value$cores)
da <- list()
if ("ANCOMBC2" %in% methods && requireNamespace("ANCOMBC", quietly = TRUE)) {
  say("ANCOM-BC2 ...")
  meta_ab <- cbind(meta_da, .sample_id = rownames(meta_da))
  ab <- ANCOMBC::ancombc2(
    data = counts_da, taxa_are_rows = TRUE, meta_data = meta_ab,
    fix_formula = deparse1(da_f[[2]]), rand_formula = if (nzchar(random_terms)) random_terms else NULL,
    p_adj_method = "BH", prv_cut = 0, lib_cut = 0, group = NULL, struc_zero = FALSE,
    neg_lb = FALSE, alpha = num("da_q"), pseudo_sens = TRUE, global = FALSE,
    n_cl = cores, verbose = FALSE)$res
  da$ANCOMBC2 <- do.call(rbind, lapply(coefs, function(cf) {
    q <- ab[[paste0("q_", cf)]]
    ss <- ab[[paste0("passed_ss_", cf)]]
    if (!is.null(ss)) q[!ss] <- 1   # require the pseudo-count sensitivity analysis to pass
    data.frame(feature = ab$taxon, coefficient = cf, effect = ab[[paste0("lfc_", cf)]], q = q)
  }))
}
if ("MaAsLin3" %in% methods && requireNamespace("maaslin3", quietly = TRUE)) {
  say("MaAsLin3 ...")
  mm3 <- meta_da
  mm3$log_depth <- log(colSums(counts_da)[rownames(mm3)])
  f3 <- paste("~", deparse1(da_f[[2]]), "+ log_depth",
              if (nzchar(random_terms)) paste("+", random_terms) else "")
  fit <- maaslin3::maaslin3(
    input_data = as.data.frame(t(counts_da)), input_metadata = mm3, output = out("maaslin3"),
    formula = f3, normalization = "TSS", transform = "LOG", augment = TRUE, standardize = FALSE,
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
    plot_summary_plot = FALSE, plot_associations = FALSE, verbosity = "ERROR", cores = cores)
  r <- fit$fit_data_abundance$results
  r <- r[r$name %in% coefs & is.na(r$error), ]
  da$MaAsLin3 <- data.frame(feature = r$feature, coefficient = r$name, effect = r$coef,
                            q = r$qval_individual)
}
if ("ALDEx2" %in% methods && requireNamespace("ALDEx2", quietly = TRUE)) {
  say("ALDEx2 (supporting; no random effects) ...")
  clr <- ALDEx2::aldex.clr(counts_da, mm, mc.samples = 128, denom = "all", verbose = FALSE)
  g <- ALDEx2::aldex.glm(clr, verbose = FALSE)
  da$ALDEx2 <- do.call(rbind, lapply(coefs, function(cf) {
    e <- g[[paste0(cf, ":Est")]]; q <- g[[paste0(cf, ":pval.padj")]]
    if (is.null(e) || is.null(q)) return(NULL)
    data.frame(feature = rownames(g), coefficient = cf, effect = e, q = q)
  }))
}
primary <- intersect(c("ANCOMBC2", "MaAsLin3"), names(da))
if (!length(primary)) stop("neither ANCOMBC nor maaslin3 is installed; no primary DA method")
if (length(primary) < 2) say("WARNING: only %s available; the consensus rests on one method", primary)
da_all <- do.call(rbind, lapply(names(da), function(m) cbind(method = m, da[[m]])))
da_all$role <- ifelse(da_all$method %in% primary, "primary", "supporting")
write_tsv(da_all, "da_all_methods.tsv")

qmax <- num("da_q")
sig <- da_all[!is.na(da_all$q) & da_all$q <= qmax & !is.na(da_all$effect), ]
key <- function(x) paste(x$feature, x$coefficient, sep = "\r")
cons <- NULL
for (k in unique(key(sig[sig$method %in% primary, ]))) {
  s <- sig[key(sig) == k, ]
  p <- s[s$method %in% primary, ]
  if (nrow(p) < length(primary) || length(unique(sign(p$effect))) != 1) next
  sup <- s[s$method == "ALDEx2", ]
  cons <- rbind(cons, data.frame(
    coefficient = p$coefficient[1], feature = p$feature[1], direction = if (p$effect[1] > 0) "up" else "down",
    ANCOMBC2_lfc = if ("ANCOMBC2" %in% p$method) p$effect[p$method == "ANCOMBC2"] else NA,
    ANCOMBC2_q = if ("ANCOMBC2" %in% p$method) p$q[p$method == "ANCOMBC2"] else NA,
    MaAsLin3_coef = if ("MaAsLin3" %in% p$method) p$effect[p$method == "MaAsLin3"] else NA,
    MaAsLin3_q = if ("MaAsLin3" %in% p$method) p$q[p$method == "MaAsLin3"] else NA,
    ALDEx2_agrees = if (!"ALDEx2" %in% names(da)) NA else
      nrow(sup) > 0 && sign(sup$effect[1]) == sign(p$effect[1])))
}
if (is.null(cons)) cons <- data.frame(coefficient = character(0), feature = character(0),
                                     direction = character(0))
cons <- cons[order(match(cons$coefficient, coefs), cons$feature), , drop = FALSE]
write_tsv(cons, "da_consensus.tsv")
say("DA consensus (%s, q <= %.2f, same direction): %d taxon-coefficient pairs",
    paste(primary, collapse = " + "), qmax, nrow(cons))

# ---- 6. Sensitivity depth ----------------------------------------------------------------

if (!is.na(sens_depth)) {
  say("Sensitivity depth %d: alpha diversity and Jaccard PERMANOVA again", sens_depth)
  invisible(fit_alpha(rarefied_alpha(sens_depth, iter), "_sensitivity"))
  if ("jaccard" %in% dist_names) invisible(permanova("jaccard", sens_depth, "_sensitivity"))
}

# ---- Settings used -------------------------------------------------------------------------

used <- data.frame(setting = names(value), value = unlist(value), source = unlist(source_of[names(value)]))
used <- rbind(used, data.frame(setting = c("random_effects", "main_depth_used", "sensitivity_depth_used"),
                               value = c(random_terms, main_depth, if (is.na(sens_depth)) "" else sens_depth),
                               source = "derived"))
write_tsv(used, "design_used.tsv")
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
say("Done. Results in %s", value$out_dir)
