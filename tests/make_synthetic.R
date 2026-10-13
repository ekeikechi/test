# Synthetic stand-in for pipeline/10_export output: 30 samples, 2 groups,
# 300 ASVs in 60 genera; 6 genera are truly enriched in "treated".
set.seed(1)
out <- commandArgs(TRUE)[1]
n_s <- 30; n_g <- 60; asv_per_g <- 5; n_a <- n_g * asv_per_g
samples <- sprintf("%02d-S%02d", 1:n_s, 1:n_s)          # leading digits test make.names quirk
grp <- rep(c("control", "treated"), each = n_s / 2)
site <- rep(c("A", "B"), times = n_s / 2)
genus_of <- rep(seq_len(n_g), each = asv_per_g)
base <- rlnorm(n_g, 3, 1.5)
lib <- round(runif(n_s, 4000, 20000))
da_genera <- 1:6
counts <- sapply(seq_len(n_s), function(j) {
  mu_g <- base * ifelse(seq_len(n_g) %in% da_genera & grp[j] == "treated", 6, 1)
  mu <- mu_g[genus_of] * runif(n_a, 0.5, 1.5)
  p <- mu / sum(mu)
  rnbinom(n_a, size = 2, mu = p * lib[j])
})
asv <- sprintf("asv%03d%s", seq_len(n_a), substr(sapply(seq_len(n_a), function(i) paste(sample(letters, 6), collapse = "")), 1, 6))
dimnames(counts) <- list(asv, samples)
phyla <- c("p__Bacillota_A", "p__Bacteroidota", "p__Pseudomonadota", "p__Actinomycetota")
tax <- sapply(seq_len(n_a), function(i) {
  g <- genus_of[i]; ph <- phyla[(g %% 4) + 1]
  lin <- c("d__Bacteria", ph, sprintf("c__C%d", g %% 8), sprintf("o__O%d", g %% 12),
           sprintf("f__F%d", g %% 20), sprintf("g__Genus%02d", g))
  if (g > 55) lin <- lin[1:5]                         # some genus-unassigned taxa
  if (i %% 3 == 0 && g <= 55) lin <- c(lin, sprintf("s__Species%03d", i))
  paste(lin, collapse = "; ")
})
writeLines(c("# Constructed from biom file",
             paste(c("#OTU ID", samples), collapse = "\t"),
             apply(cbind(asv, counts), 1, paste, collapse = "\t")),
           file.path(out, "feature-table.tsv"))
write.table(data.frame(`Feature ID` = asv, Taxon = tax, Confidence = 0.95, check.names = FALSE),
            file.path(out, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
md <- data.frame(`sample-id` = samples, treatment = grp, site = site,
                 pH = round(rnorm(n_s, 7, 0.4), 2), check.names = FALSE)
md_lines <- c(paste(names(md), collapse = "\t"), "#q2:types\tcategorical\tcategorical\tnumeric",
              apply(md, 1, paste, collapse = "\t"))
writeLines(md_lines, file.path(out, "metadata.tsv"))
tr <- ape::rtree(n_a + 40, rooted = TRUE, tip.label = c(asv, sprintf("ref%02d", 1:40)))
ape::write.tree(tr, file.path(out, "tree.nwk"))
cat("truth: enriched genera =", sprintf("Genus%02d", da_genera), "\n")
