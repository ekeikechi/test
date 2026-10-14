#!/usr/bin/env Rscript
# =============================================================================
# PR2 trait annotation of ASVs (stage 11a; runs on the stage 10 export).
#
# Usage:
#   Rscript R/pr2_traits.R [EXPORT_DIR] [OUT_DIR] [PR2_FILE]
#     EXPORT_DIR  stage 10 output (feature-table.tsv, taxonomy.tsv, metadata.tsv,
#                 analysis_settings.tsv); default results/10_export
#     OUT_DIR     default results/11_traits
#     PR2_FILE    the full PR2 database (pr2.qs from the pr2database R package, or an
#                 .rds of the same list); default: download it once into OUT_DIR
#
# Needs: R with the 'qs' package (conda-forge r-qs) to read pr2.qs; ggplot2 for plots.
#
# What it does
#   1. Traits of taxa (PR2 taxonomy table, one row per species), with the overrides in
#      R/pr2_trait_overrides.tsv applied first:
#        ecological_function (Sommeria-Klein et al. 2021), mixoplankton (Mitra et al.
#        2023), HAB_species (WoRMS/IOC-UNESCO HAB list), WoRMS habitat flags, WoRMS and
#        GBIF ids.
#   2. Each ASV takes the values of the PR2 species under its deepest assigned taxon
#      (matched on the full lineage). A categorical trait is assigned when at least
#      trait_min_agreement (default 0.9) of those species share one value; otherwise it
#      is "mixed". Every value is reported with the rank used and the number of species.
#   3. photo_category, from ecological_function + mixoplankton + class:
#        plant                     Embryophyceae (land plants), always reported apart
#        acquired_phototrophy      mixoplankton pSNCM, GNCM or eSNCM
#        phototroph                ecological_function "phototrophs", or mixoplankton CM
#        dinoflagellate_unresolved ecological_function "dinoflagellates" without CM
#        non_photosynthetic        phagotrophs, parasites, metazoans, osmotrophs
#        mixed_or_unknown          no agreement, or no PR2 trait data for the taxon
#        unassigned                ASV not classified (no taxon)
#   4. Evidence columns (ev_*): what the PR2 reference sequences of the taxon record about
#      their source (biotic relationship, hosts, biome/material, habitat keywords,
#      share from cultures, geography). Descriptive, sparse and biased toward well-sampled
#      environments; not traits.
#   5. Cross-references: consensus SILVA and EukRibo/UniEuk taxonomy of the same sequences.
#   6. habitat_flag: with sample_habitat set (marine, brackish, freshwater, terrestrial),
#      a taxon with WoRMS habitat data for >= 3 species but none in the sample habitat is
#      flagged for checking (marine and brackish count as each other). Not a filter.
#      WoRMS habitat columns read "k of n": k of the n species with WoRMS habitat data
#      have a record in that habitat (a species can have several).
#   Consensus SILVA/UniEuk strings are majority-rule: the deepest prefix shared by more
#   than half of the taxon's reference sequences.
#
# Outputs (OUT_DIR): asv_traits.tsv, table_by_photo_category.tsv,
#   table_by_ecological_function.tsv (+ _relative), trait_summary.tsv,
#   traits_settings.tsv, photo_category_by_sample.pdf, phototroph_share_by_group.pdf,
#   export_phototroph/ (stage 10-style export of phototroph ASVs only, for
#   R/downstream_analysis.R).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg <- list(
  export_dir = if (length(args) >= 1) args[[1]] else "results/10_export",
  out_dir    = if (length(args) >= 2) args[[2]] else "results/11_traits",
  pr2_file   = if (length(args) >= 3) args[[3]] else NULL,
  pr2_url    = "https://storage.googleapis.com/pr2database-data/pr2database/data-qs/pr2.qs",
  overrides  = NULL,                 # default: pr2_trait_overrides.tsv next to this script
  ranks      = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order",
                 "Family", "Genus", "Species"),
  min_agreement  = 0.9,              # TRAIT_MIN_AGREEMENT
  sample_habitat = "",               # SAMPLE_HABITAT
  pr2_version    = "",               # PR2_VERSION of the classifier (for a version check)
  group          = NULL,             # GROUP_COL
  habitat_min_species = 3
)
say <- function(...) message(sprintf(...))
PR2_LEVELS <- c("domain", "supergroup", "division", "subdivision", "class", "order",
                "family", "genus", "species")

# ---- settings handed over by stage 10 -----------------------------------------
st_path <- file.path(cfg$export_dir, "analysis_settings.tsv")
if (file.exists(st_path)) {
  st <- read.delim(st_path, colClasses = "character", quote = "", na.strings = character(0))
  st <- setNames(st$value, st$setting)
  get <- function(k) if (k %in% names(st) && nzchar(st[[k]])) st[[k]] else NULL
  if (!is.null(v <- get("taxonomy_ranks")))      cfg$ranks <- strsplit(v, "\\s+")[[1]]
  if (!is.null(v <- get("trait_min_agreement"))) cfg$min_agreement <- as.numeric(v)
  if (!is.null(v <- get("sample_habitat")))      cfg$sample_habitat <- v
  if (!is.null(v <- get("pr2_version")))         cfg$pr2_version <- v
  if (!is.null(v <- get("group_col")))           cfg$group <- v
}
if (length(cfg$ranks) != 9)
  stop(sprintf("PR2 traits need the 9 PR2 ranks; analysis_settings.tsv lists %d: %s",
               length(cfg$ranks), paste(cfg$ranks, collapse = " ")))
if (nzchar(cfg$sample_habitat) &&
    !cfg$sample_habitat %in% c("marine", "brackish", "freshwater", "terrestrial"))
  stop("sample_habitat must be marine, brackish, freshwater, terrestrial or empty")
if (is.null(cfg$overrides)) {
  here <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  cfg$overrides <- file.path(if (is.na(here)) "R" else dirname(here), "pr2_trait_overrides.tsv")
}
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(...) file.path(cfg$out_dir, ...)

# ---- PR2 database -------------------------------------------------------------
if (is.null(cfg$pr2_file)) {
  cfg$pr2_file <- out("pr2.qs")
  if (!file.exists(cfg$pr2_file)) {
    say("Downloading the PR2 database (about 160 MB) from %s", cfg$pr2_url)
    options(timeout = max(1800, getOption("timeout")))
    download.file(cfg$pr2_url, cfg$pr2_file, mode = "wb", quiet = TRUE)
  }
}
pr2 <- if (grepl("\\.rds$", cfg$pr2_file, ignore.case = TRUE)) readRDS(cfg$pr2_file) else {
  if (!requireNamespace("qs", quietly = TRUE))
    stop("Reading pr2.qs needs the R package 'qs' (e.g. conda install -c conda-forge r-qs).")
  qs::qread(cfg$pr2_file)
}
pr2_version <- if (!is.null(pr2$version)) as.character(pr2$version) else "unknown"
say("PR2 database version %s (%s)", pr2_version, cfg$pr2_file)
if (nzchar(cfg$pr2_version) && cfg$pr2_version != pr2_version)
  warning(sprintf(paste("The classifier was built from PR2 %s but the trait data are PR2 %s.",
                        "Taxa renamed between versions will get no traits."),
                  cfg$pr2_version, pr2_version), call. = FALSE)

tax <- as.data.frame(pr2$taxonomy, stringsAsFactors = FALSE)
tax <- tax[tax$domain %in% "Eukaryota", , drop = FALSE]
main <- as.data.frame(pr2$main, stringsAsFactors = FALSE)
main <- main[main$gene %in% "18S_rRNA" & main$domain %in% "Eukaryota" &
               (is.na(main$organelle) | main$organelle %in% "nucleus"), , drop = FALSE]
for (col in c("mixoplankton", "HAB_species", "ecological_function", "worms_id", "gbif_id",
              "worms_marine", "worms_brackish", "worms_freshwater", "worms_terrestrial"))
  if (!col %in% names(tax)) tax[[col]] <- NA
for (col in c("eukref_biotic_relationship", "eukref_specific_host", "gb_host",
              "eukref_env_biome", "eukref_env_material", "gb_isolation_source",
              "pr2_sample_type", "pr2_ocean", "pr2_continent", "silva_taxonomy",
              "eukribo_UniEuk_taxonomy_string"))
  if (!col %in% names(main)) main[[col]] <- NA_character_
say("PR2 traits: %d eukaryote taxa (species), %d nuclear 18S reference sequences",
    nrow(tax), nrow(main))

# ---- overrides ----------------------------------------------------------------
ov <- read.delim(cfg$overrides, colClasses = "character", quote = "")
tax$override <- ""
for (i in seq_len(nrow(ov))) {
  lvl <- tolower(ov$rank[i])
  if (!lvl %in% PR2_LEVELS) stop(sprintf("override %d: unknown rank '%s'", i, ov$rank[i]))
  if (!ov$field[i] %in% names(tax)) stop(sprintf("override %d: unknown field '%s'", i, ov$field[i]))
  hit <- which(tax[[lvl]] == ov$taxon[i])
  if (!length(hit)) { say("Override not used (taxon not in PR2 %s): %s %s", pr2_version, ov$rank[i], ov$taxon[i]); next }
  tax[[ov$field[i]]][hit] <- ov$value[i]
  tax$override[hit] <- paste0(tax$override[hit], ifelse(nzchar(tax$override[hit]), "; ", ""),
                              sprintf("%s %s: %s=%s", ov$rank[i], ov$taxon[i], ov$field[i], ov$value[i]))
}
say("Overrides applied: %d of %d rules (%s)", sum(ov$taxon %in% unlist(tax[PR2_LEVELS])),
    nrow(ov), cfg$overrides)

# ---- lineage keys -------------------------------------------------------------
lineage_keys <- function(m) {          # m: character matrix, one column per level
  keys <- vector("list", ncol(m)); acc <- m[, 1]
  keys[[1]] <- acc
  for (k in seq_len(ncol(m))[-1]) { acc <- paste(acc, m[, k], sep = ";"); keys[[k]] <- acc }
  keys
}
tax_keys  <- lineage_keys(as.matrix(tax[PR2_LEVELS]))
main_keys <- lineage_keys(as.matrix(main[PR2_LEVELS]))
tax_index  <- lapply(tax_keys,  function(k) split(seq_along(k), k))
main_index <- lapply(main_keys, function(k) split(seq_along(k), k))

# ---- ASV taxonomy ---------------------------------------------------------------
tx <- read.delim(file.path(cfg$export_dir, "taxonomy.tsv"), check.names = FALSE,
                 comment.char = "", stringsAsFactors = FALSE, quote = "")
labels <- lapply(strsplit(tx$Taxon, ";\\s*"), function(x) {
  x <- trimws(x); x <- x[nzchar(x)]
  x <- sub("^[a-z]+__", "", x)
  cut <- which(!nzchar(x) | x %in% c("Unassigned", "unassigned"))
  if (length(cut)) x <- x[seq_len(cut[1] - 1)]
  x
})
depth <- lengths(labels)
if (any(depth > 9)) stop("taxonomy.tsv is deeper than the 9 PR2 ranks")

# ---- summaries ------------------------------------------------------------------
top_value <- function(x, min_agree) {          # categorical: assigned value + agreement
  x <- x[!is.na(x) & nzchar(x) & x != "unknown"]
  if (!length(x)) return(list(value = NA_character_, agreement = NA_real_, n = 0L))
  t <- sort(table(x), decreasing = TRUE); share <- t[[1]] / length(x)
  list(value = if (share >= min_agree) names(t)[1] else "mixed", agreement = round(share, 3),
       n = length(x))
}
fmt_counts <- function(x, k = 4) {
  x <- trimws(x[!is.na(x)]); x <- x[nzchar(x)]
  if (!length(x)) return("")
  t <- sort(table(x), decreasing = TRUE)[seq_len(min(k, length(unique(x))))]
  paste(sprintf("%s (%d)", names(t), as.integer(t)), collapse = "; ")
}
consensus <- function(x, sep) {   # majority-rule: deepest prefix shared by > 50% of sequences
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(c("", ""))
  parts <- strsplit(x, sep, fixed = TRUE); best <- c("", "")
  for (j in seq_len(max(lengths(parts)))) {
    pre <- vapply(parts, function(p) if (length(p) >= j) paste(p[1:j], collapse = sep) else NA_character_, "")
    t <- sort(table(pre[!is.na(pre)]), decreasing = TRUE)
    if (!length(t) || t[[1]] / length(x) <= 0.5) break
    best <- c(names(t)[1], sprintf("%.2f", t[[1]] / length(x)))
  }
  best
}
habitat_patterns <- c(            # keyword classes for ev_habitat_keywords (case-insensitive)
  marine     = "marine|seawater|sea water|\\bsea\\b|ocean|coastal|offshore|pelagic|\\bbay\\b|gulf|fjord",
  brackish   = "brackish|estuar|lagoon|salt ?marsh",
  freshwater = "fresh ?water|\\blake|river|\\bpond|stream|reservoir|creek|wetland|\\bbog\\b|\\bpeat",
  soil       = "\\bsoil|rhizosphere|leaf litter|compost|grassland|forest floor",
  sediment   = "sediment|\\bmud|tidal flat|intertidal|\\bsand|benthic|beach",
  host       = "fece|faece|fecal|faecal|\\bgut\\b|intestin|rumen|ruminal|blood|tissue|stool|infect|lesion|\\bskin")
habitat_keywords <- function(rows) {
  txt <- tolower(paste(main$gb_isolation_source[rows], main$eukref_env_biome[rows],
                       main$eukref_env_material[rows], sep = " | "))
  host <- nzchar(paste0(ifelse(is.na(main$gb_host[rows]), "", main$gb_host[rows]),
                        ifelse(is.na(main$eukref_specific_host[rows]), "", main$eukref_specific_host[rows])))
  n <- vapply(names(habitat_patterns), function(h) {
    hit <- grepl(habitat_patterns[[h]], txt, perl = TRUE)
    if (h == "host") hit <- hit | host
    sum(hit)
  }, numeric(1))
  n <- n[n > 0]
  if (!length(n)) "" else paste(sprintf("%s (%d)", names(n), as.integer(n))[order(-n)], collapse = "; ")
}
hab_count <- function(v) sum(v %in% 1)   # species with a WoRMS record in that habitat

annotate <- function(lab) {
  d <- length(lab)
  res <- list(assigned_rank = if (d) cfg$ranks[d] else "", assigned_taxon = if (d) lab[d] else "")
  if (!d) return(c(res, photo_category = "unassigned"))
  key <- paste(lab, collapse = ";")
  # deepest level (<= d) of this lineage present in PR2's taxonomy table
  lvl <- d; rows <- tax_index[[d]][[key]]
  while (is.null(rows) && lvl > 1) { lvl <- lvl - 1; rows <- tax_index[[lvl]][[paste(lab[1:lvl], collapse = ";")]] }
  if (is.null(rows)) return(c(res, trait_rank = "", n_pr2_species = 0L, photo_category = "mixed_or_unknown"))
  t <- tax[rows, , drop = FALSE]
  eco <- top_value(t$ecological_function, cfg$min_agreement)
  mixo <- top_value(ifelse(is.na(t$mixoplankton), "not_listed", t$mixoplankton), cfg$min_agreement)
  hab_sp <- t$species[t$HAB_species %in% "HAB species"]
  w <- t[!is.na(t$worms_marine) | !is.na(t$worms_brackish) | !is.na(t$worms_freshwater) |
           !is.na(t$worms_terrestrial), , drop = FALSE]
  counts <- c(marine = hab_count(w$worms_marine), brackish = hab_count(w$worms_brackish),
              freshwater = hab_count(w$worms_freshwater), terrestrial = hab_count(w$worms_terrestrial))
  of_n <- function(h) sprintf("%d of %d", counts[[h]], nrow(w))
  flag <- ""
  if (nzchar(cfg$sample_habitat) && nrow(w) >= cfg$habitat_min_species) {
    target <- switch(cfg$sample_habitat, marine = c("marine", "brackish"),
                     brackish = c("brackish", "marine"), cfg$sample_habitat)
    if (sum(counts[target]) == 0)
      flag <- sprintf("no %s WoRMS record among %d species with habitat data; check assignment",
                      cfg$sample_habitat, nrow(w))
  }
  plant <- d >= 5 && lab[5] == "Embryophyceae"
  photo <- if (plant) "plant"
    else if (mixo$value %in% c("pSNCM", "GNCM", "eSNCM")) "acquired_phototrophy"
    else if (eco$value %in% "phototrophs" || mixo$value %in% "CM") "phototroph"
    else if (eco$value %in% "dinoflagellates") "dinoflagellate_unresolved"
    else if (eco$value %in% c("phagotrophs", "parasites", "metazoans", "osmotrophs")) "non_photosynthetic"
    else "mixed_or_unknown"
  m_rows <- main_index[[lvl]][[paste(lab[1:lvl], collapse = ";")]]
  if (is.null(m_rows)) m_rows <- integer(0)
  st_type <- main$pr2_sample_type[m_rows]; st_type <- st_type[!is.na(st_type)]
  silva <- consensus(main$silva_taxonomy[m_rows], ";")
  unieuk <- consensus(main$eukribo_UniEuk_taxonomy_string[m_rows], "|")
  sp_level <- lvl == 9 && nrow(t) == 1
  c(res, list(
    trait_rank = cfg$ranks[lvl], n_pr2_species = nrow(t),
    ecological_function = eco$value, ecological_function_agreement = eco$agreement,
    ecological_function_n = eco$n,
    mixoplankton = mixo$value, mixoplankton_agreement = mixo$agreement,
    photo_category = photo,
    HAB = if (sp_level) (if (length(hab_sp)) "HAB species" else "no")
          else if (length(hab_sp)) sprintf("taxon contains %d HAB species", length(hab_sp)) else "no",
    HAB_species = paste(head(hab_sp, 5), collapse = "; "),
    worms_n_species = nrow(w),
    worms_marine = of_n("marine"), worms_brackish = of_n("brackish"),
    worms_freshwater = of_n("freshwater"), worms_terrestrial = of_n("terrestrial"),
    habitat_flag = flag,
    worms_id = if (sp_level && !is.na(t$worms_id)) as.character(t$worms_id) else "",
    gbif_id = if (sp_level && !is.na(t$gbif_id)) as.character(t$gbif_id) else "",
    overrides_applied = paste(unique(t$override[nzchar(t$override)]), collapse = "; "),
    ev_n_ref_seqs = length(m_rows),
    ev_biotic_relationship = fmt_counts(main$eukref_biotic_relationship[m_rows]),
    ev_hosts = fmt_counts(c(main$eukref_specific_host[m_rows], main$gb_host[m_rows]), 3),
    ev_env_biome_material = fmt_counts(c(main$eukref_env_biome[m_rows], main$eukref_env_material[m_rows]), 4),
    ev_habitat_keywords = if (length(m_rows)) habitat_keywords(m_rows) else "",
    ev_culture_share = if (length(st_type)) round(mean(st_type %in% c("culture", "isolate")), 3) else NA,
    ev_geography = fmt_counts(c(main$pr2_ocean[m_rows], main$pr2_continent[m_rows]), 4),
    silva_consensus = silva[1], silva_consensus_share = silva[2],
    unieuk_consensus = unieuk[1], unieuk_consensus_share = unieuk[2]))
}

ukeys <- vapply(labels, paste, character(1), collapse = ";")
uniq <- !duplicated(ukeys)
say("Annotating %d ASVs (%d distinct lineages) ...", nrow(tx), sum(uniq))
ann_u <- lapply(labels[uniq], annotate)
cols <- unique(unlist(lapply(ann_u, names)))
to_df <- function(a) as.data.frame(lapply(setNames(cols, cols), function(cn) {
  v <- lapply(a, function(x) if (is.null(x[[cn]])) NA else x[[cn]])
  unlist(lapply(v, function(z) if (length(z) == 1) z else NA))
}), stringsAsFactors = FALSE)
ann <- to_df(ann_u)[match(ukeys, ukeys[uniq]), , drop = FALSE]
traits <- cbind(feature_id = tx[[1]], taxon = tx$Taxon, ann, stringsAsFactors = FALSE)

# ---- feature table ----------------------------------------------------------------
ft_path <- file.path(cfg$export_dir, "feature-table.tsv")
ft_lines <- readLines(ft_path)
ft <- read.delim(ft_path, skip = 1, check.names = FALSE, comment.char = "", row.names = 1,
                 stringsAsFactors = FALSE)
ft <- as.matrix(ft)
traits$total_count <- as.numeric(rowSums(ft))[match(traits$feature_id, rownames(ft))]
write.table(traits, out("asv_traits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

collapse_by <- function(field) {
  g <- traits[[field]][match(rownames(ft), traits$feature_id)]
  g[is.na(g) | !nzchar(g)] <- "no_data"
  agg <- rowsum(ft, g)
  write.table(data.frame(category = rownames(agg), agg, check.names = FALSE),
              out(sprintf("table_by_%s.tsv", field)), sep = "\t", quote = FALSE, row.names = FALSE)
  rel <- sweep(agg, 2, pmax(colSums(agg), 1), "/")
  write.table(data.frame(category = rownames(rel), round(rel, 5), check.names = FALSE),
              out(sprintf("table_by_%s_relative.tsv", field)), sep = "\t", quote = FALSE, row.names = FALSE)
  rel
}
rel_photo <- collapse_by("photo_category")
invisible(collapse_by("ecological_function"))

summ <- do.call(rbind, lapply(split(traits, traits$photo_category), function(d)
  data.frame(photo_category = d$photo_category[1], n_asvs = nrow(d),
             reads = sum(d$total_count, na.rm = TRUE))))
summ$read_share <- round(summ$reads / sum(summ$reads), 4)
summ <- summ[order(-summ$reads), ]
write.table(summ, out("trait_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Photosynthesis categories (all samples pooled):")
message(paste(sprintf("  %-26s %6d ASVs  %6.1f%% of reads", summ$photo_category, summ$n_asvs,
                      100 * summ$read_share), collapse = "\n"))
nflag <- sum(nzchar(traits$habitat_flag) & !is.na(traits$habitat_flag))
if (nzchar(cfg$sample_habitat))
  say("Habitat check (%s): %d ASVs flagged (see habitat_flag in asv_traits.tsv)", cfg$sample_habitat, nflag)

write.table(data.frame(
  setting = c("pr2_trait_version", "pr2_file", "pr2_url", "classifier_pr2_version",
              "overrides_file", "overrides_rules", "trait_min_agreement", "sample_habitat",
              "habitat_min_species", "date"),
  value = c(pr2_version, cfg$pr2_file, cfg$pr2_url, cfg$pr2_version, cfg$overrides, nrow(ov),
            cfg$min_agreement, cfg$sample_habitat, cfg$habitat_min_species, format(Sys.time()))),
  out("traits_settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- phototroph-only export (input for R/downstream_analysis.R) ---------------------
sub_dir <- out("export_phototroph")
dir.create(sub_dir, showWarnings = FALSE)
keep_ids <- traits$feature_id[traits$photo_category == "phototroph"]
keep <- c(1, 2, 2 + which(sub("\t.*", "", ft_lines[-(1:2)]) %in% keep_ids))
writeLines(ft_lines[keep], file.path(sub_dir, "feature-table.tsv"))
write.table(tx[tx[[1]] %in% keep_ids, , drop = FALSE], file.path(sub_dir, "taxonomy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (f in c("metadata.tsv", "analysis_settings.tsv"))
  if (file.exists(file.path(cfg$export_dir, f)))
    file.copy(file.path(cfg$export_dir, f), file.path(sub_dir, f), overwrite = TRUE)
say("Phototroph-only export: %d ASVs in %s", length(keep_ids), sub_dir)

# ---- plots ---------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  long <- data.frame(sample = rep(colnames(rel_photo), each = nrow(rel_photo)),
                     category = rep(rownames(rel_photo), ncol(rel_photo)),
                     share = as.vector(rel_photo))
  md_path <- file.path(cfg$export_dir, "metadata.tsv")
  md <- if (file.exists(md_path)) read.delim(md_path, check.names = FALSE, colClasses = "character",
                                             comment.char = "") else NULL
  has_group <- !is.null(md) && !is.null(cfg$group) && cfg$group %in% names(md)
  if (has_group) long$group <- md[[cfg$group]][match(long$sample, md[[1]])]
  p <- ggplot2::ggplot(long, ggplot2::aes(sample, share, fill = category)) +
    ggplot2::geom_col(width = 1) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 4)) +
    ggplot2::labs(x = NULL, y = "Share of reads", fill = "Photosynthesis category")
  if (has_group) p <- p + ggplot2::facet_grid(~group, scales = "free_x", space = "free_x")
  ggplot2::ggsave(out("photo_category_by_sample.pdf"), p, width = 12, height = 5)
  if (has_group) {
    ph <- long[long$category == "phototroph", ]
    q <- ggplot2::ggplot(ph, ggplot2::aes(group, share)) + ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.15, size = 1) + ggplot2::theme_bw() +
      ggplot2::labs(x = cfg$group, y = "Phototroph share of reads")
    ggplot2::ggsave(out("phototroph_share_by_group.pdf"), q, width = 6, height = 4)
  }
}
say("Done: %s", cfg$out_dir)
