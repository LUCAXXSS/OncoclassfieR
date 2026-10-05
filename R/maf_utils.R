# MAF harmonisation, oncoplot and TMB helpers -----------------------------

#' Common column names shared by several MAF data frames
#'
#' Returns the columns that are present in all MAF data frames of a list.
#' Useful before merging MAF data coming from different cohorts.
#'
#' @param maf_list List of MAF objects or MAF-like data frames.
#'
#' @return Character vector of common column names.
#'
#' @export
#' @examples
#' common <- maf_common_columns(list(example_maf, example_maf))
#' head(common)
maf_common_columns <- function(maf_list) {
  cols <- lapply(maf_list, function(m) {
    colnames(.maf_to_dataframe(m))
  })
  Reduce(intersect, cols)
}

#' Harmonise a MAF data frame to a common column set
#'
#' Restricts a MAF data frame to the columns shared across cohorts, tags each
#' row with a `cohort_name`, and optionally keeps only a gene panel. This
#' mirrors the harmonisation step used to merge WES / panel cohorts before
#' Jaccard clustering.
#'
#' @param maf MAF object or MAF-like data frame.
#' @param common_columns Character vector of columns to keep (typically the
#'   output of [maf_common_columns()]).
#' @param gene_use Character vector of genes to keep; `NULL` keeps all genes.
#' @param cohort_name Character tag stored in the new `cohort_name` column.
#'
#' @return A harmonised data frame.
#'
#' @export
#' @examples
#' common <- maf_common_columns(list(example_maf))
#' shaped <- shape_maf(example_maf, common, gene_use = NULL, cohort_name = "WES")
#' head(colnames(shaped))
shape_maf <- function(maf, common_columns, gene_use = NULL, cohort_name = "cohort") {
  maf <- .maf_to_dataframe(maf)
  maf_small <- maf[, colnames(maf) %in% common_columns, drop = FALSE]

  if (!is.null(cohort_name)) {
    maf_small$cohort_name <- cohort_name
  }

  if (!is.null(gene_use)) {
    maf_small <- maf_small[maf_small$Hugo_Symbol %in% gene_use, , drop = FALSE]
  }

  maf_small
}

#' Oncoplot for one sample cluster
#'
#' Subsets a MAF object to the samples of a single cluster and draws an
#' oncoplot of the `top` most frequent genes in that cluster.
#'
#' @param maf A `MAF` object.
#' @param cluster_assignment Result list of [jaccard_cluster()] or a data
#'   frame with `Tumor_Sample_Barcode` and `cluster` columns.
#' @param cluster_id Cluster to plot (as character or numeric).
#' @param top Integer. Number of genes shown by the oncoplot.
#' @param colors Named character vector of variant-class colours (passed to
#'   [maftools::oncoplot()]).
#' @param ... Further arguments passed to [maftools::oncoplot()].
#'
#' @return Invisibly returns the filtered `MAF` object; the oncoplot is drawn
#'   on the current graphics device.
#'
#' @export
#' @examples
#' \donttest{
#' res <- jaccard_cluster(example_maf, n_clusters = 3, plot = FALSE)
#' maf <- maftools::read.maf(example_maf)
#' plot_cluster_oncoplot(maf, res, cluster_id = 1, top = 5)
#' }
plot_cluster_oncoplot <- function(maf,
                                  cluster_assignment,
                                  cluster_id,
                                  top = 10,
                                  colors = NULL,
                                  ...) {
  if (is.list(cluster_assignment) && !is.data.frame(cluster_assignment) &&
      "sample_clusters" %in% names(cluster_assignment)) {
    cluster_assignment <- cluster_assignment$sample_clusters
  }

  selected <- cluster_assignment$Tumor_Sample_Barcode[
    as.character(cluster_assignment$cluster) == as.character(cluster_id)
  ]

  if (length(selected) == 0) {
    stop("No samples found for cluster `", cluster_id, "`.", call. = FALSE)
  }

  maf_filtered <- maftools::subsetMaf(maf, tsb = selected)

  title <- paste0("Cluster ", cluster_id, " (n=", length(selected), ")")
  maftools::oncoplot(maf_filtered, top = top, titleText = title,
                     colors = colors, ...)

  invisible(maf_filtered)
}

#' Boxplot of TMB across sample clusters
#'
#' Computes TMB with [maftools::tmb()], merges it with the cluster
#' assignment, and draws a boxplot with jittered points plus pairwise
#' Wilcoxon tests (Bonferroni adjusted).
#'
#' @param maf A `MAF` object.
#' @param cluster_assignment Result list of [jaccard_cluster()] or a data
#'   frame with `Tumor_Sample_Barcode` and `cluster` columns.
#' @param comparisons Optional list of cluster pairs to annotate, e.g.
#'   `list(c("1", "3"))`. If `NULL`, all pairs are tested.
#'
#' @return A list with the `ggplot` object and the pairwise Wilcoxon test
#'   table.
#'
#' @export
#' @examples
#' \donttest{
#' res <- jaccard_cluster(example_maf, n_clusters = 3, plot = FALSE)
#' maf <- maftools::read.maf(example_maf)
#' out <- plot_tmb_by_cluster(maf, res)
#' out$wilcox
#' }
plot_tmb_by_cluster <- function(maf, cluster_assignment, comparisons = NULL) {
  if (!methods::is(maf, "MAF")) {
    stop("`maf` must be a MAF object.", call. = FALSE)
  }
  if (is.list(cluster_assignment) && !is.data.frame(cluster_assignment) &&
      "sample_clusters" %in% names(cluster_assignment)) {
    cluster_assignment <- cluster_assignment$sample_clusters
  }

  sample_tmb <- maftools::tmb(maf)
  info <- merge(cluster_assignment, sample_tmb,
                by = "Tumor_Sample_Barcode", all.x = TRUE)

  info$cluster <- as.character(info$cluster)

  p <- ggplot2::ggplot(
    info,
    ggplot2::aes(x = .data$cluster, y = log(.data$total_perMB + 1, 2),
                 fill = .data$cluster)
  ) +
    ggplot2::geom_boxplot(outliers = FALSE) +
    ggplot2::geom_jitter(color = "gray30", width = 0.1, alpha = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Cluster", y = "log2(TMB + 1)", fill = "Cluster") +
    ggsci::scale_fill_npg()

  if (!is.null(comparisons)) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons, method = "wilcox"
    )
  }

  wilcox_res <- NULL
  if (length(unique(info$cluster)) >= 2) {
    pw <- stats::pairwise.wilcox.test(
      log(info$total_perMB + 1, 2), info$cluster,
      p.adjust.method = "bonferroni",
      exact = FALSE
    )
    wilcox_res <- as.data.frame(pw$p.value)
  }

  list(plot = p, wilcox = wilcox_res)
}
