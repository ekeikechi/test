# Synthetic stand-in for a split-plot field experiment, laid out like the Grande Anse design:
#   3 rounds x 2 nutrient levels x 3 replicate plots = 18 plots (whole plots; each plot is
#   sampled in one round), 3 exclosure subplots per plot, 4 sediment depths per subplot
#   = 216 samples, 40 genera x 3 ASVs, 9-rank PR2-style taxonomy.
# Known effects (multipliers on expected abundance):
#   Genus01-03  x6 in nutrient "enriched"          (whole-plot stratum)
#   Genus04-06  x6 in exclosure "full_exclosure"   (subplot stratum)
#   Genus07-09  x1, 0.4, 0.15, 0.05 at 0/1/4/7 cm  (sample stratum)
# Every genus also gets a random plot effect (sd 0.6 on the log scale) and a smaller random
# subplot effect (sd 0.3), so a model that ignores the plots treats correlated samples as
# independent. Two samples (GA1UC7S04, GA2FB4S08) are shallow (600 reads), below
# MIN_SAMPLE_DEPTH = 1000.
#
# Usage: Rscript tests/make_synthetic_splitplot.R OUT_DIR
set.seed(42)
out <- commandArgs(TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

plots <- expand.grid(rep = 1:3, nutrient = c("F", "U"), round = 1:3, stringsAsFactors = FALSE)
plots$plot <- sprintf("S%02d", seq_len(nrow(plots)))
d <- merge(plots, expand.grid(exclosure = c("B", "C", "F"), depth_cm = c(0, 1, 4, 7),
                              stringsAsFactors = FALSE))
d <- d[order(d$round, d$plot, d$exclosure, d$depth_cm), ]
d$sample <- sprintf("GA%d%s%s%d%s", d$round, d$nutrient, d$exclosure, d$depth_cm, d$plot)
n_s <- nrow(d); n_g <- 40; asv_per_g <- 3; n_a <- n_g * asv_per_g
genus_of <- rep(seq_len(n_g), each = asv_per_g)

base <- rlnorm(n_g, 3, 1.2)
plot_eff <- matrix(rnorm(length(unique(d$plot)) * n_g, 0, 0.6), ncol = n_g,
                   dimnames = list(unique(d$plot), NULL))
sub_id <- paste(d$plot, d$exclosure, sep = "_")
sub_eff <- matrix(rnorm(length(unique(sub_id)) * n_g, 0, 0.3), ncol = n_g,
                  dimnames = list(unique(sub_id), NULL))
depth_mult <- c(`0` = 1, `1` = 0.4, `4` = 0.15, `7` = 0.05)
lib <- round(runif(n_s, 3000, 15000))
lib[d$sample %in% c("GA1UC7S04", "GA2FB4S08")] <- 600

counts <- sapply(seq_len(n_s), function(j) {
  mu_g <- base * exp(plot_eff[d$plot[j], ] + sub_eff[sub_id[j], ])
  if (d$nutrient[j] == "F") mu_g[1:3] <- mu_g[1:3] * 6
  if (d$exclosure[j] == "F") mu_g[4:6] <- mu_g[4:6] * 6
  mu_g[7:9] <- mu_g[7:9] * depth_mult[[as.character(d$depth_cm[j])]]
  mu <- mu_g[genus_of] * runif(n_a, 0.5, 1.5)
  rnbinom(n_a, size = 3, mu = mu / sum(mu) * lib[j])
})
asv <- sprintf("asv%03d", seq_len(n_a))
dimnames(counts) <- list(asv, d$sample)

tax <- vapply(seq_len(n_a), function(i) {
  g <- genus_of[i]
  paste0(paste(c("d__Eukaryota", "sgr__TSAR", sprintf("dv__Division%d", g %% 3),
                 sprintf("dvs__Subdivision%d", g %% 4), sprintf("c__C%d", g %% 6),
                 sprintf("o__O%d", g %% 8), sprintf("f__F%d", g %% 12),
                 sprintf("g__Genus%02d", g)), collapse = "; "), ";")
}, character(1))

writeLines(c("# Constructed from biom file",
             paste(c("#OTU ID", d$sample), collapse = "\t"),
             apply(cbind(asv, counts), 1, paste, collapse = "\t")),
           file.path(out, "feature-table.tsv"))
write.table(data.frame(`Feature ID` = asv, Taxon = tax, Confidence = 0.95, check.names = FALSE),
            file.path(out, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lab <- c(B = "bird_exclosure", C = "open_control", F = "full_exclosure")
nut <- c(F = "enriched", U = "unenriched_control")
md <- data.frame(`sample-id` = d$sample, round = paste0("R", d$round),
                 round_label = c("early_July", "late_July", "late_August")[d$round],
                 nutrient = nut[d$nutrient], exclosure = lab[d$exclosure],
                 depth_cm = d$depth_cm, depth_layer = paste0(d$depth_cm, "cm"),
                 plot = sub("^S", "P", d$plot),
                 subplot = paste(sub("^S", "P", d$plot), d$exclosure, sep = "_"),
                 treatment = paste(nut[d$nutrient], lab[d$exclosure], sep = "__"),
                 check.names = FALSE)
writeLines(c(paste(names(md), collapse = "\t"),
             paste(c("#q2:types", rep("categorical", 4), "numeric", rep("categorical", 4)),
                   collapse = "\t"),
             apply(md, 1, paste, collapse = "\t")),
           file.path(out, "metadata.tsv"))
writeLines(c("setting\tvalue", "preset\t18S_PR2",
             "taxonomy_ranks\tDomain Supergroup Division Subdivision Class Order Family Genus Species",
             "collapse_level\t8", "overview_level\t3", "sampling_depth\tauto",
             "group_col\ttreatment",
             "da_reference\tnutrient::unenriched_control exclosure::open_control",
             "random_seed\t711", "min_sample_depth\t1000",
             "depth_sensitivity_groups\tround depth_layer", "depth_sensitivity_min_keep\t0.8"),
           file.path(out, "analysis_settings.tsv"))
cat("truth: nutrient -> Genus01-03, exclosure (full) -> Genus04-06, depth -> Genus07-09;",
    "shallow samples GA1UC7S04 GA2FB4S08\n")
