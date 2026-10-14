# Synthetic PR2 database (pr2.rds: list(main, taxonomy, version)) and a stage 10-style
# export with 14 ASVs whose trait annotations are known. Used by tests/traits_check.sh.
# Usage: Rscript tests/make_synthetic_pr2_traits.R OUTDIR [CLASSIFIER_PR2_VERSION]
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[[1]]; clf_version <- if (length(args) >= 2) args[[2]] else "5.1.1"
dir.create(file.path(outdir, "export"), recursive = TRUE, showWarnings = FALSE)
lv <- c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")

sp <- function(path, species, eco, mixo = NA, hab = NA, wm = NA, wb = NA, wf = NA, wt = NA,
               worms = NA, gbif = NA) {
  p <- strsplit(path, ";")[[1]]
  d <- as.data.frame(as.list(setNames(c(p, species), lv)), stringsAsFactors = FALSE)
  cbind(d, mixoplankton = mixo, HAB_species = hab, ecological_function = eco,
        worms_id = worms, worms_marine = wm, worms_brackish = wb, worms_freshwater = wf,
        worms_terrestrial = wt, gbif_id = gbif, n_seq_SSU = 1, stringsAsFactors = FALSE)
}
NAV <- "Eukaryota;TSAR;Stramenopiles;Gyrista;Bacillariophyceae;Naviculales;Naviculaceae;Navicula"
GYR <- "Eukaryota;TSAR;Stramenopiles;Gyrista;Bacillariophyceae;Naviculales;Pleurosigmataceae;Gyrosigma"
CHR <- "Eukaryota;TSAR;Stramenopiles;Gyrista;Chrysophyceae;Chromulinales"
PRO <- "Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Prorocentrales;Prorocentraceae;Prorocentrum"
PER <- "Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Peridiniales;Protoperidiniaceae;Protoperidinium"
MES <- "Eukaryota;TSAR;Alveolata;Ciliophora;Litostomatea;Cyclotrichiida;Mesodiniidae;Mesodinium"
MIC <- "Eukaryota;Archaeplastida;Chlorophyta;Chlorophyta_X;Trebouxiophyceae;Chlorellales;Chlorellaceae;Micractinium"
SPA <- "Eukaryota;Archaeplastida;Streptophyta;Streptophyta_X;Embryophyceae;Poales;Poaceae;Spartina"
GYM <- "Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae"
CIL <- "Eukaryota;TSAR;Alveolata;Ciliophora;Spirotrichea;Strombidiida"
tax <- rbind(
  sp(NAV, "Navicula_a", "phototrophs", wm = 1, wb = 0, wf = 0, wt = 0),
  sp(NAV, "Navicula_b", "phototrophs", wm = 1, wb = 1, wf = 0, wt = 0),
  sp(NAV, "Navicula_c", "phototrophs"),
  sp(GYR, "Gyrosigma_fasciola", "phototrophs", wm = 0, wb = 0, wf = 1, wt = 0),
  sp(paste0(CHR, ";Paraphysomonadaceae;Paraphysomonas"), "Paraphysomonas_vestita", "phototrophs"),
  sp(paste0(CHR, ";Paraphysomonadaceae;Paraphysomonas"), "Paraphysomonas_foraminifera", "phototrophs"),
  sp(paste0(CHR, ";Ochromonadaceae;Ochromonas"), "Ochromonas_danica", "phototrophs", mixo = "CM"),
  sp(paste0(CHR, ";Ochromonadaceae;Ochromonas"), "Ochromonas_sp.", "phototrophs"),
  sp(paste0(CHR, ";Chromulinaceae;Spumella"), "Spumella_elongata", "phagotrophs"),
  sp(paste0(CHR, ";Chromulinaceae;Spumella"), "Spumella_sp.", "phagotrophs"),
  sp(PRO, "Prorocentrum_minimum", "dinoflagellates", mixo = "CM", hab = "HAB species",
     wm = 1, worms = 12345, gbif = 678),
  sp(PRO, "Prorocentrum_lima", "dinoflagellates", mixo = "CM", hab = "HAB species", wm = 1),
  sp(PER, "Protoperidinium_a", "dinoflagellates"),
  sp(PER, "Protoperidinium_b", "dinoflagellates"),
  sp(MES, "Mesodinium_rubrum", "phagotrophs", mixo = "pSNCM", wm = 1),
  sp(MIC, "Micractinium_a", "phototrophs", wm = 0, wb = 0, wf = 1, wt = 0),
  sp(MIC, "Micractinium_b", "phototrophs", wm = 0, wb = 0, wf = 1, wt = 0),
  sp(MIC, "Micractinium_c", "phototrophs", wm = 0, wb = 0, wf = 1, wt = 1),
  sp(SPA, "Spartina_alterniflora", "phototrophs", wm = 1),
  sp(paste0(GYM, ";Gymnodinium"), "Gymnodinium_a", "dinoflagellates"),
  sp(paste0(GYM, ";Amphidinium"), "Amphidinium_a", "dinoflagellates"),
  sp(paste0(GYM, ";Unlisteddinium"), "Unlisteddinium_a", "dinoflagellates"),
  sp(paste0(CIL, ";Strombidiidae;Strombidium"), "Strombidium_a", "phagotrophs"),
  sp(paste0(CIL, ";Unlistedidae;Unlistedia"), "Unlistedia_a", "phagotrophs"))

# reference sequences: one per species, three for Navicula_a with source records
main <- tax[, lv]; main$gene <- "18S_rRNA"; main$organelle <- "nucleus"
for (cn in c("eukref_biotic_relationship", "eukref_specific_host", "gb_host", "eukref_env_biome",
             "eukref_env_material", "gb_isolation_source", "pr2_sample_type", "pr2_ocean",
             "pr2_continent", "silva_taxonomy", "eukribo_UniEuk_taxonomy_string"))
  main[[cn]] <- NA_character_
nav <- main[rep(which(main$species == "Navicula_a"), 2), ]
main <- rbind(main, nav)
i <- which(main$species == "Navicula_a")
main$gb_isolation_source[i] <- c("intertidal mudflat sediment, Bay of Fundy", "culture", NA)
main$gb_host[i] <- c(NA, NA, "Caretta caretta")
main$pr2_sample_type[i] <- c("environmental", "culture", "culture")
main$eukref_biotic_relationship[i] <- c("free-living", NA, NA)
main$pr2_ocean[i] <- c("Atlantic Ocean", "Atlantic Ocean", NA)
main$silva_taxonomy[i] <- c("Eukaryota;SAR;Diatomea;Navicula;Navicula a", "Eukaryota;SAR;Diatomea;Navicula;Navicula b",
                            "Eukaryota;SAR;Diatomea;Navicula;uncultured")
saveRDS(list(main = main, taxonomy = tax, version = "5.1.1"), file.path(outdir, "pr2.rds"))

pre <- c("d__", "sgr__", "dv__", "dvs__", "c__", "o__", "f__", "g__", "s__")
lin <- function(path, species = NULL) {
  p <- c(strsplit(path, ";")[[1]], species)
  paste(paste0(pre[seq_along(p)], p), collapse = "; ")
}
asv <- c(
  asv01 = lin(NAV),                                   # phototroph, genus
  asv02 = lin(paste0(CHR, ";Paraphysomonadaceae;Paraphysomonas"), "Paraphysomonas_vestita"),  # override
  asv03 = paste0(lin("Eukaryota;TSAR;Stramenopiles;Gyrista;Chrysophyceae"), "; o__; f__"),    # mixed class
  asv04 = lin(SPA, "Spartina_alterniflora"),          # plant
  asv05 = lin(PRO),                                   # CM dinoflagellates -> phototroph; HAB genus
  asv06 = lin(PER),                                   # dinoflagellate_unresolved
  asv07 = lin(MES, "Mesodinium_rubrum"),              # acquired_phototrophy
  asv08 = lin(MIC),                                   # habitat flag (3 freshwater-only species)
  asv09 = lin(GYR, "Gyrosigma_fasciola"),             # 1 species: no flag
  asv10 = "Unassigned",
  asv11 = lin(NAV, "Navicula_newsp"),                 # species not in PR2 -> genus
  asv12 = "d__Eukaryota; sgr__TSAR; dv__Unknowndivision",   # -> supergroup, mixed
  asv13 = lin(PRO, "Prorocentrum_minimum"),           # HAB species, WoRMS/GBIF ids
  asv14 = lin(NAV, "Navicula_a"),                     # evidence columns
  asv15 = lin(paste0(GYM, ";Gymnodinium")),           # Adl Table 2.1 L M -> phototroph
  asv16 = lin(paste0(GYM, ";Amphidinium")),           # Adl: species-level variation -> unresolved
  asv17 = lin(paste0(GYM, ";Unlisteddinium")),        # not in Adl table -> unresolved
  asv18 = lin(paste0(CIL, ";Strombidiidae;Strombidium")),  # Adl Table 2.2 genus row: O
  asv19 = lin(paste0(CIL, ";Unlistedidae;Unlistedia")))    # ciliate default B
write.table(data.frame(`Feature ID` = names(asv), Taxon = asv, Confidence = 0.95, check.names = FALSE),
            file.path(outdir, "export", "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
smp <- sprintf("S%d", 1:6)
m <- matrix(10L * seq_along(asv), nrow = length(asv), ncol = 6, dimnames = list(names(asv), smp))
con <- file(file.path(outdir, "export", "feature-table.tsv"), "w")
writeLines(c("# Constructed from biom file", paste(c("#OTU ID", smp), collapse = "\t")), con)
write.table(m, con, sep = "\t", quote = FALSE, col.names = FALSE); close(con)
write.table(data.frame(`sample-id` = smp, zone = rep(c("upper", "lower"), each = 3), check.names = FALSE),
            file.path(outdir, "export", "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("setting\tvalue",
             "taxonomy_ranks\tDomain Supergroup Division Subdivision Class Order Family Genus Species",
             "group_col\tzone", "sample_habitat\tmarine", "trait_min_agreement\t0.9",
             paste0("pr2_version\t", clf_version)),
           file.path(outdir, "export", "analysis_settings.tsv"))
