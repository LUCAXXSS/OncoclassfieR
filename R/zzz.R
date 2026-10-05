# Package-level documentation -------------------------------------------

#' OncoclassfieR: Jaccard Similarity Based Sample Classification for Cancer MAF Data
#'
#' Classifies tumour samples into molecular subgroups based on the Jaccard
#' similarity of their mutation profiles. The package is built around the
#' sample-clustering workflow used in the LCBM / HSYY panel analysis:
#' \itemize{
#'   \item Pre-process MAF objects (\code{\link{preprocess_maf_top_genes}}).
#'   \item Compute pairwise Jaccard similarity and cluster samples
#'     (\code{\link{jaccard_similarity_matrix}}, \code{\link{jaccard_cluster}}).
#'   \item Re-draw the similarity heatmap in a given subgroup order
#'     (\code{\link{plot_jaccard_heatmap}}).
#'   \item Rank genes per cluster and visualise them
#'     (\code{\link{cluster_top_genes}}, \code{\link{plot_cluster_piechart}}).
#'   \item Run downstream analyses: oncoplots per cluster
#'     (\code{\link{plot_cluster_oncoplot}}), TMB comparison
#'     (\code{\link{plot_tmb_by_cluster}}) and survival analysis
#'     (\code{\link{fit_km}}, \code{\link{plot_km_curves}},
#'     \code{\link{fit_cox}}).
#' }
#'
#' @docType package
#' @name OncoclassfieR
#' @keywords internal
"_PACKAGE"

# Quiet R CMD check notes about tidy-eval / NSE used in plotting helpers.
utils::globalVariables(c(
  "cluster", "count", "gene", "gene_label", "percentage", "total_count",
  "isTreated", "Value", "isTreat", "cohort", "Percentage", "total_perMB",
  "Tumor_Sample_Barcode", ".data"
))
