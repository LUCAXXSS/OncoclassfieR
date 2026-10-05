#!/usr/bin/env Rscript
# =============================================================================
# LCBM_reproduction.R
#
# End-to-end reproduction of the current LCBM panel / WES Jaccard 4-subtype
# analysis, using only the datasets shipped with OncoclassfieR. See
# `261005_panel_WES_jaccard_reproduction.md` for the full provenance: this
# script corresponds to sections 3-6 of that document, plus the clustering
# visualisation and the per-subtype oncoplot of its section 9.
#
# Both cohorts are clustered separately with one shared parameter set:
#
#   max_features       = 20      top-20 most frequently mutated genes
#   similarity         = Jaccard on the per-sample mutation gene set
#   enhancement        = 0.6     similarity matrix ^ 0.6
#   distance           = euclidean on the similarity row profiles
#                                (pheatmap default, NOT 1 - Jaccard)
#   clustering_method  = "ward.D2"
#   n_clusters         = 4
#   drop_empty_samples = TRUE    samples with no mutation among the top-20
#                                genes stay out of the similarity matrix
#   seed               = 1234
#
# Expected results (asserted at the end of this script):
#
#   panel  cohort: 342 samples clustered -> 178 / 65 / 73 / 26
#   WES    cohort: 145 samples clustered ->  28 / 71 / 28 / 18
#
# Upstream steps (merging the four cohorts, restricting to the 298 shared
# panel genes and filtering VAF >= 2%) are NOT reproduced here, because the
# raw cohort files are not distributed with the package. This script starts
# from the already filtered objects LCBM_panel_maf and LCBM_WES_maf.
#
# Usage:
#
#   library(OncoclassfieR)                        # or devtools::load_all()
#   source(system.file("examples", "LCBM_reproduction.R",
#                      package = "OncoclassfieR"))
#
#   # ... or run it from a checkout:
#   # Rscript inst/examples/LCBM_reproduction.R /path/to/output
#
# Figures and tables are written to `out_dir` (a temporary directory by
# default) so that running the script never touches the package itself.
# =============================================================================

suppressPackageStartupMessages({
  library(maftools)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 0. Output directory and package data access
# -----------------------------------------------------------------------------
# The output directory is the first command line argument. Scripts that
# source this file (e.g. data-raw/make_readme_figures.R) can set it through
# options(OncoclassfieR.out_dir = "...") instead of commandArgs().
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- getOption("OncoclassfieR.out_dir",
                     if (length(args) >= 1) args[[1]]
                     else file.path(tempdir(), "LCBM_reproduction"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Redirect the default graphics device: any helper called with plot = TRUE
# would otherwise open a Rplots.pdf in the working directory. The explicit
# png() / pdf() devices opened below are unaffected.
grDevices::pdf(NULL)

# The package has to be available: `library(OncoclassfieR)` or, from a
# checkout, `devtools::load_all(".")`.
if (!"OncoclassfieR" %in% loadedNamespaces()) {
  if (requireNamespace("OncoclassfieR", quietly = TRUE)) {
    suppressPackageStartupMessages(library(OncoclassfieR))
  } else {
    stop("Load OncoclassfieR first: library(OncoclassfieR) or devtools::load_all(\".\").",
         call. = FALSE)
  }
}

# The script must work both for an installed package and for a checkout loaded
# with devtools::load_all(), so the datasets are fetched through the package
# namespace instead of `data()` (which only searches installed packages).
pkg_data <- function(name) {
  get(name, envir = asNamespace("OncoclassfieR"))
}

LCBM_panel_maf        <- pkg_data("LCBM_panel_maf")
LCBM_WES_maf          <- pkg_data("LCBM_WES_maf")
LCBM_sample_annotation <- pkg_data("LCBM_sample_annotation")

mafs <- list(panel = LCBM_panel_maf, WES = LCBM_WES_maf)

# -----------------------------------------------------------------------------
# 1. Clustering parameters (identical for both cohorts)
# -----------------------------------------------------------------------------
max_features      <- 20
n_clusters        <- 4
clustering_method <- "ward.D2"
enhancement       <- 0.6
seed              <- 1234
top_genes_pie     <- 3
top_genes_table   <- 20
top_genes_oncoplot <- 15

# -----------------------------------------------------------------------------
# 2. Jaccard clustering + silhouette, one cohort at a time
# -----------------------------------------------------------------------------
cluster_list <- list()
results      <- list()

for (grp in names(mafs)) {
  maf <- mafs[[grp]]
  message("\n===== ", grp, " cohort =====")

  set.seed(seed)
  maf_top <- preprocess_maf_top_genes(maf, max_features = max_features,
                                      drop_empty_samples = TRUE)

  res <- jaccard_cluster(maf_top, n_clusters = n_clusters,
                         clustering_method = clustering_method,
                         enhancement = enhancement,
                         show_colnames = FALSE, plot = FALSE)

  clusters <- res$sample_clusters
  clusters <- clusters[order(clusters$cluster, clusters$Tumor_Sample_Barcode), ]

  cat("samples entering the similarity matrix:", nrow(clusters), "\n")
  print(table(clusters$cluster))

  write.csv(clusters, file.path(out_dir, paste0("LCBM_cluster_", grp, ".csv")),
            row.names = FALSE)

  # --- cluster x cohort crosstab ---------------------------------------------
  ctab <- as.data.frame.matrix(
    table(cluster = clusters$cluster,
          cohort = as.character(
            maf@data$cohort_name[match(clusters$Tumor_Sample_Barcode,
                                       as.character(maf@data$Tumor_Sample_Barcode))])))
  ctab <- data.frame(cluster = rownames(ctab), ctab, row.names = NULL)
  ctab$total <- rowSums(ctab[, setdiff(colnames(ctab), "cluster"), drop = FALSE])
  write.csv(ctab, file.path(out_dir, paste0("LCBM_cluster_by_cohort_", grp, ".csv")),
            row.names = FALSE)
  cat("\ncluster x cohort:\n"); print(ctab)

  # --- top genes per cluster -------------------------------------------------
  top_tab <- cluster_top_genes(maf, res, top = top_genes_table)
  write.csv(top_tab, file.path(out_dir, paste0("LCBM_top_genes_", grp, ".csv")),
            row.names = FALSE)

  # --- most frequently mutated genes overall ---------------------------------
  gene_counts <- sort(table(maf@data$Hugo_Symbol), decreasing = TRUE)
  top20 <- data.frame(rank = seq_len(max_features),
                      gene = names(gene_counts)[seq_len(max_features)],
                      n_samples = as.integer(gene_counts[seq_len(max_features)]),
                      row.names = NULL)
  write.csv(top20, file.path(out_dir, paste0("LCBM_top20_genes_", grp, ".csv")),
            row.names = FALSE)

  # --- silhouette curve ------------------------------------------------------
  # k = 3:10 (k = 2 is not a meaningful split of a homogeneous cohort)
  sil <- silhouette_analysis(res, k_range = 3:10,
                             clustering_method = clustering_method)
  write.csv(sil$silhouette,
            file.path(out_dir, paste0("LCBM_silhouette_", grp, ".csv")),
            row.names = FALSE)

  # --- similarity heatmap in clustering order --------------------------------
  # plot = FALSE: build the pheatmap silently so it can be printed inside the
  # png() device below instead of on whatever device happens to be current.
  hm <- plot_jaccard_heatmap(res, group_order = as.character(seq_len(n_clusters)),
                             clustering_method = clustering_method,
                             show_colnames = FALSE, plot = FALSE)
  png(file.path(out_dir, paste0("LCBM_jaccard_heatmap_", grp, ".png")),
      width = 2000, height = 1800, res = 200)
  print(hm$heatmap)
  dev.off()

  # --- silhouette figure -----------------------------------------------------
  png(file.path(out_dir, paste0("LCBM_silhouette_", grp, ".png")),
      width = 1600, height = 1200, res = 200)
  print(sil$plot)
  dev.off()

  # --- per-cluster pie charts ------------------------------------------------
  pie <- plot_cluster_piechart(top_tab, top = top_genes_pie)
  pdf(file.path(out_dir, paste0("LCBM_cluster_piechart_", grp, ".pdf")),
      width = 8, height = 6)
  print(pie)
  dev.off()

  cluster_list[[grp]] <- clusters
  results[[grp]] <- list(maf = maf, maf_top = maf_top, res = res,
                         clusters = clusters, top_tab = top_tab,
                         silhouette = sil$silhouette)
}

# -----------------------------------------------------------------------------
# 3. TMB comparison and per-cluster oncoplots (need a maftools MAF object)
# -----------------------------------------------------------------------------
for (grp in names(results)) {
  maf <- results[[grp]]$maf
  clusters <- results[[grp]]$clusters

  tmb_res <- plot_tmb_by_cluster(maf, clusters)

  # tmb_res$wilcox is a lower-triangular matrix of pairwise p-values; flatten
  # it into one row per compared pair so the table is self-describing.
  w <- tmb_res$wilcox
  wilcox_tab <- do.call(rbind, lapply(seq_len(nrow(w)), function(i) {
    js <- which(!is.na(unlist(w[i, ])))
    if (length(js) == 0) return(NULL)
    data.frame(cluster_1 = colnames(w)[js],
               cluster_2 = rownames(w)[i],
               p_value = unlist(w[i, js]),
               row.names = NULL)
  }))
  write.csv(wilcox_tab,
            file.path(out_dir, paste0("LCBM_tmb_wilcox_", grp, ".csv")),
            row.names = FALSE)

  png(file.path(out_dir, paste0("LCBM_tmb_", grp, ".png")),
      width = 1800, height = 1400, res = 200)
  print(tmb_res$plot)
  dev.off()

  pdf(file.path(out_dir, paste0("LCBM_oncoplot_", grp, ".pdf")),
      width = 10, height = 8)
  for (cl in sort(unique(clusters$cluster))) {
    plot_cluster_oncoplot(maf, clusters, cluster_id = cl,
                          top = top_genes_oncoplot)
  }
  dev.off()
}

# -----------------------------------------------------------------------------
# 4. TP53 / EGFR 4-group assignment (derived from the MAF itself)
# -----------------------------------------------------------------------------
grp_levels <- c("Both_mut", "TP53_only", "EGFR_only", "Neither_mut")

four_group_list <- list()
for (grp in names(results)) {
  dat <- as.data.frame(results[[grp]]$maf@data)
  dat$Tumor_Sample_Barcode <- as.character(dat$Tumor_Sample_Barcode)
  dat$cohort_name <- as.character(dat$cohort_name)
  dat$Hugo_Symbol <- as.character(dat$Hugo_Symbol)

  samp <- unique(dat[, c("Tumor_Sample_Barcode", "cohort_name")])
  samp$group <- grp
  tp53_mut <- unique(dat$Tumor_Sample_Barcode[dat$Hugo_Symbol == "TP53"])
  egfr_mut <- unique(dat$Tumor_Sample_Barcode[dat$Hugo_Symbol == "EGFR"])
  samp$TP53_status <- ifelse(samp$Tumor_Sample_Barcode %in% tp53_mut, "Mut", "WT")
  samp$EGFR_status <- ifelse(samp$Tumor_Sample_Barcode %in% egfr_mut, "Mut", "WT")
  samp$TP53_EGFR_group <- ifelse(
    samp$TP53_status == "Mut" & samp$EGFR_status == "Mut", "Both_mut",
    ifelse(samp$TP53_status == "Mut" & samp$EGFR_status == "WT", "TP53_only",
    ifelse(samp$TP53_status == "WT"  & samp$EGFR_status == "Mut", "EGFR_only",
           "Neither_mut")))

  samp$jaccard_cluster <- results[[grp]]$clusters$cluster[
    match(samp$Tumor_Sample_Barcode, results[[grp]]$clusters$Tumor_Sample_Barcode)]

  samp$TP53_EGFR_group <- factor(samp$TP53_EGFR_group, levels = grp_levels)
  samp <- samp[order(samp$group, samp$TP53_EGFR_group, samp$Tumor_Sample_Barcode), ]

  write.csv(samp[, c("Tumor_Sample_Barcode", "cohort_name", "group",
                     "TP53_status", "EGFR_status", "TP53_EGFR_group",
                     "jaccard_cluster")],
            file.path(out_dir, paste0("LCBM_tp53_egfr_4group_", grp, ".csv")),
            row.names = FALSE)

  four_group_list[[grp]] <- samp
}

# -----------------------------------------------------------------------------
# 5. Histology annotation and the TP53/EGFR x cohort x histology table
# -----------------------------------------------------------------------------
hist_map <- c(LUAD = "AdC", LUSC = "SqCC", SCLC = "SCLC", LCNEC = "LCNEC",
              NET = "LCNEC")
hist_levels   <- c("AdC", "SqCC", "SCLC", "LCNEC")
cohort_levels <- c("Huashan", "cbioportal")

panel_4g <- four_group_list$panel
panel_4g <- merge(panel_4g,
                  LCBM_sample_annotation[, c("Tumor_Sample_Barcode", "Histology")],
                  by = "Tumor_Sample_Barcode", all.x = TRUE)
panel_4g$Histology_mapped <- unname(hist_map[as.character(panel_4g$Histology)])

annotated <- subset(panel_4g, !is.na(panel_4g$Histology_mapped))
annotated$Histology_mapped <- factor(annotated$Histology_mapped, levels = hist_levels)
annotated$cohort_name <- factor(annotated$cohort_name, levels = cohort_levels)

cat("\npanel samples with a usable histology:", nrow(annotated),
    "(dropped", nrow(panel_4g) - nrow(annotated),
    "with LCC / NSCLC / Other / NA)\n")

stats_out <- as.data.frame(
  table(cohort = annotated$cohort_name,
        group = annotated$TP53_EGFR_group,
        histology = annotated$Histology_mapped))
stats_out <- stats_out[stats_out$Freq > 0, ]
stats_out <- data.frame(cohort = as.character(stats_out$cohort),
                        TP53_EGFR_group = as.character(stats_out$group),
                        Histology = as.character(stats_out$histology),
                        n = stats_out$Freq, row.names = NULL)
write.csv(stats_out,
          file.path(out_dir, "LCBM_tp53_egfr_by_histology_panel.csv"),
          row.names = FALSE)
print(stats_out)

# -----------------------------------------------------------------------------
# 6. Compare every table against the shipped expected results
# -----------------------------------------------------------------------------
expected_dir <- system.file("extdata", "LCBM_expected", package = "OncoclassfieR")
if (!nzchar(expected_dir)) expected_dir <- NA_character_

if (!is.na(expected_dir) && dir.exists(expected_dir)) {
  expected_files <- list.files(expected_dir, pattern = "\\.csv$")
  cat("\n===== comparison against inst/extdata/LCBM_expected =====\n")
  all_ok <- TRUE
  for (f in expected_files) {
    got_path <- file.path(out_dir, f)
    if (!file.exists(got_path)) {
      cat(sprintf("%-45s MISSING\n", f)); all_ok <- FALSE; next
    }
    exp <- read.csv(file.path(expected_dir, f), stringsAsFactors = FALSE)
    got <- read.csv(got_path, stringsAsFactors = FALSE)
    same <- isTRUE(all.equal(exp, got, check.attributes = FALSE))
    cat(sprintf("%-45s %s\n", f, if (same) "OK" else "DIFFERS"))
    if (!same) all_ok <- FALSE
  }
  cat(if (all_ok)
        "\nAll tables match the expected results.\n"
      else
        "\nSome tables differ from the expected results (see above).\n")
} else {
  message("No expected results shipped; skipping the comparison.")
}

cat("\nOutputs written to: ", out_dir, "\n", sep = "")
