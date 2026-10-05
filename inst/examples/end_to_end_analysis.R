# =============================================================================
# OncoclassfieR: end-to-end example analysis
# -----------------------------------------------------------------------------
# Runs the complete LCBM / HSYY-style workflow on the bundled example data:
#   1. Pre-process the MAF (keep the most frequently mutated genes)
#   2. Compute the pairwise Jaccard similarity matrix
#   3. Cluster samples by hierarchical clustering
#   4. Rank the top genes per cluster and draw pie charts
#   5. Merge the discovered clusters with survival data
#   6. Kaplan-Meier and Cox survival analysis
#
# Usage (run from the package root):
#   Rscript inst/examples/end_to_end_analysis.R
#
# All plots and tables are written to ./output_end_to_end.
# =============================================================================

## 0. Setup ----------------------------------------------------------------
# Run this script from the package root. When the checkout is present it is
# loaded from source with devtools::load_all(); otherwise the installed
# package is used.

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(OncoclassfieR))
}

out_dir <- "output_end_to_end"
dir.create(out_dir, showWarnings = FALSE)

n_clusters <- 4L        # number of molecular subgroups to look for
max_features <- 20L     # keep the 20 most frequently mutated genes

## 1. Input data ------------------------------------------------------------
# Simulated MAF-like mutation data (40 samples, 4 latent mutation programmes)
# and matching overall-survival data shipped with the package.
data(example_maf)
data(example_survival)

cat("example_maf:", nrow(example_maf), "mutations in",
    length(unique(example_maf$Tumor_Sample_Barcode)), "samples\n")

   ## 2. Pre-processing ---------------------------------------------------------
# Rare genes would dominate the pairwise similarities, so reduce the MAF to
# the genes mutated most often across the cohort.
maf_top <- preprocess_maf_top_genes(example_maf, max_features = max_features)
cat("Genes kept after pre-processing:",
    length(unique(maf_top$Hugo_Symbol)), "\n")

## 3. Jaccard similarity matrix ----------------------------------------------
sim <- jaccard_similarity_matrix(maf_top, enhancement = 0.6)
cat("Jaccard similarity matrix:", nrow(sim), "x", ncol(sim), "\n")
write.csv(round(sim, 3), file.path(out_dir, "jaccard_similarity.csv"))

## 4. Sample clustering -------------------------------------------------------
# Clusters the samples and draws the similarity heatmap.
png(file.path(out_dir, "jaccard_heatmap.png"),
    width = 8, height = 6, units = "in", res = 300)
res <- jaccard_cluster(maf_top, n_clusters = n_clusters,
                       clustering_method = "ward.D2")
dev.off()

cat("Cluster sizes:\n")
print(table(res$sample_clusters$cluster))
write.csv(res$sample_clusters,
          file.path(out_dir, "sample_clusters.csv"), row.names = FALSE)

## 4b. Heatmap in a hand-picked subgroup order --------------------------------
# Figures usually have to keep a fixed subgroup order (to line up with the
# oncoplots, survival curves, ...). plot_jaccard_heatmap() re-draws the
# similarity matrix in exactly that order, leaving res untouched.
cluster_order <- c(3, 1, 4, 2)

png(file.path(out_dir, "jaccard_heatmap_ordered.png"),
    width = 8, height = 6, units = "in", res = 300)
hm <- plot_jaccard_heatmap(res, group_order = cluster_order,
                           show_colnames = FALSE)
dev.off()

stopifnot(identical(hm$group_order, as.character(cluster_order)))
write.csv(hm$group, file.path(out_dir, "sample_clusters_ordered.csv"),
          row.names = FALSE)

## 5. Top genes per cluster ---------------------------------------------------
cluster_genes <- cluster_top_genes(maf_top, res, top = 3)
print(cluster_genes)

p_pie <- plot_cluster_piechart(cluster_genes, top = 3)
ggplot2::ggsave(file.path(out_dir, "cluster_top_genes_piechart.png"),
                p_pie, width = 10, height = 5, dpi = 300)

## 6. Survival data merged with the discovered clusters -----------------------
# The bundled survival table uses the ground-truth cluster labels; merge the
# *discovered* cluster assignment so downstream analyses reflect the actual
# clustering result.
# Drop the ground-truth `cluster` column (used only to simulate the data)
# so the merge keeps the clusters *discovered* by jaccard_cluster().
surv_data <- merge(
  example_survival[, setdiff(colnames(example_survival), "cluster")],
  res$sample_clusters,
  by.x = "ID", by.y = "Tumor_Sample_Barcode", all.x = TRUE
)
stopifnot(nrow(surv_data) == nrow(example_survival))

## 7. Kaplan-Meier analysis ----------------------------------------------------
fit <- fit_km(surv_data, time = "OS", event = "status", group = "cluster")
print(fit)

km_plot <- plot_km_curves(fit, surv_data)
png(file.path(out_dir, "km_curves.png"), width = 7, height = 6,
    units = "in", res = 300)
print(km_plot)
dev.off()

## 8. Cox proportional hazards model --------------------------------------------
cox_fit <- fit_cox(surv_data, time = "OS", event = "status", group = "cluster")
print(summary(cox_fit))

## 9. Real MAF data (optional) --------------------------------------------------
# With a real MAF object the same functions work directly; uncomment and adapt:
#
#   maf <- maftools::read.maf("path/to/your.maf")
#   maf_top <- preprocess_maf_top_genes(maf, max_features = 50)
#   res <- jaccard_cluster(maf_top, n_clusters = n_clusters,
#                          clustering_method = "ward.D2")
#   plot_cluster_oncoplot(maf, res, cluster_id = 1, top = 10)  # per-cluster oncoplot
#   tmb_res <- plot_tmb_by_cluster(maf, res)                   # TMB boxplot
#   print(tmb_res$wilcox)

cat("\nAll outputs written to:", out_dir, "\n")
