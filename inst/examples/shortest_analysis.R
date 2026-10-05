# =============================================================================
# OncoclassfieR: shortest path to a real analysis
# -----------------------------------------------------------------------------
# The smallest script that runs the Jaccard 4-subtype workflow on a real
# cohort: the LCBM whole-exome cohort shipped with the package, clustered with
# the parameters of the current analysis (top-20 genes, enhancement 0.6,
# ward.D2, k = 4, empty samples dropped). It then adds a per-cluster oncoplot,
# the TMB comparison and the histology breakdown of the clusters.
#
# Usage (run from the package root):
#   Rscript inst/examples/shortest_analysis.R
#
# See LCBM_reproduction.R for the full end-to-end reproduction of both the
# panel and the WES cohorts.
# =============================================================================

## 0. Setup ----------------------------------------------------------------
if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(OncoclassfieR))
}

suppressPackageStartupMessages(library(maftools))

out_dir <- "output_shortest_analysis"
dir.create(out_dir, showWarnings = FALSE)

## 1. Input data ------------------------------------------------------------
# Real LCBM WES cohort (157 samples) and the clinical annotation table.
maf <- OncoclassfieR::LCBM_WES_maf
annotation <- OncoclassfieR::LCBM_sample_annotation

cat("LCBM WES cohort:", nrow(maf@data), "variants in",
    length(unique(maf@data$Tumor_Sample_Barcode)), "samples\n")

## 2. Clustering ------------------------------------------------------------
set.seed(1234)
maf_top <- preprocess_maf_top_genes(maf, max_features = 20,
                                    drop_empty_samples = TRUE)
res <- jaccard_cluster(maf_top, n_clusters = 4, enhancement = 0.6,
                       clustering_method = "ward.D2",
                       show_colnames = FALSE, plot = FALSE)

print(table(res$sample_clusters$cluster))   # 28 / 71 / 28 / 18

## 3. Per-cluster oncoplot and TMB ------------------------------------------
pdf(file.path(out_dir, "oncoplot_cluster1.pdf"), width = 10, height = 8)
plot_cluster_oncoplot(maf, res, cluster_id = 1, top = 20)
dev.off()

tmb_res <- plot_tmb_by_cluster(maf, res)
print(tmb_res$wilcox)

## 4. Histology of the clusters ---------------------------------------------
# The sample barcode prefixes are meaningless in the harmonised object, so the
# histology comes from the bundled annotation table instead.
maf_cluster_table <- res$sample_clusters
maf_cluster_table$Histology <- annotation$Histology[
  match(maf_cluster_table$Tumor_Sample_Barcode,
        annotation$Tumor_Sample_Barcode)]

print(table(maf_cluster_table$cluster, maf_cluster_table$Histology,
            useNA = "ifany"))

write.csv(maf_cluster_table,
          file.path(out_dir, "sample_clusters_WES.csv"), row.names = FALSE)
cat("\nAll outputs written to:", out_dir, "\n")
