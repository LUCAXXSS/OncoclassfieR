# Per-cluster gene ranking and visualisation ------------------------------

#' Top genes per sample cluster
#'
#' For each cluster returned by [jaccard_cluster()], counts how often each
#' gene is mutated among the samples of that cluster and returns the
#' `top` most frequent genes.
#'
#' @param df Data frame (or `MAF` object) with `Tumor_Sample_Barcode` and
#'   `Hugo_Symbol` columns.
#' @param cluster_assignment Either the full result list from
#'   [jaccard_cluster()], or a data frame with columns
#'   `Tumor_Sample_Barcode` and `cluster`.
#' @param top Integer. Number of top genes to return per cluster.
#'
#' @return A data frame with columns `cluster`, `gene` and `count`.
#'
#' @export
#' @examples
#' res <- jaccard_cluster(example_maf, n_clusters = 3, plot = FALSE)
#' genes <- cluster_top_genes(example_maf, res, top = 5)
#' head(genes)
cluster_top_genes <- function(df, cluster_assignment, top = 10) {
  df <- .maf_to_dataframe(df)
  df <- .dedup_sample_gene(df)

  if (is.list(cluster_assignment) && !is.data.frame(cluster_assignment) &&
      "sample_clusters" %in% names(cluster_assignment)) {
    cluster_assignment <- cluster_assignment$sample_clusters
  }
  if (!is.data.frame(cluster_assignment) ||
      !all(c("Tumor_Sample_Barcode", "cluster") %in% colnames(cluster_assignment))) {
    stop("`cluster_assignment` must be a data frame with columns ",
         "`Tumor_Sample_Barcode` and `cluster`, or the result list of ",
         "`jaccard_cluster()`.", call. = FALSE)
  }

  cluster_samples <- split(cluster_assignment$Tumor_Sample_Barcode,
                           cluster_assignment$cluster)

  cluster_genes <- list()
  for (cl in names(cluster_samples)) {
    samples <- cluster_samples[[cl]]
    df_s <- df[df$Tumor_Sample_Barcode %in% samples, , drop = FALSE]

    gene_counts <- table(df_s$Hugo_Symbol)
    gene_counts <- sort(gene_counts, decreasing = TRUE)
    n_keep <- min(top, length(gene_counts))

    cluster_genes[[cl]] <- data.frame(
      cluster = cl,
      gene = names(gene_counts)[seq_len(n_keep)],
      count = as.integer(unname(gene_counts))[seq_len(n_keep)],
      row.names = NULL
    )
  }

  do.call(rbind, cluster_genes)
}

#' Pie chart of the top genes per cluster
#'
#' Draws one pie chart per cluster showing the percentage contributed by the
#' `top` most frequent genes; the remaining mutations are pooled into an
#' "other" slice.
#'
#' @param cluster_genes Output of [cluster_top_genes()]: a data frame with
#'   columns `cluster`, `gene` and `count`.
#' @param top Integer. Number of top genes to show per cluster.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @examples
#' res <- jaccard_cluster(example_maf, n_clusters = 3, plot = FALSE)
#' genes <- cluster_top_genes(example_maf, res, top = 3)
#' p <- plot_cluster_piechart(genes, top = 3)
plot_cluster_piechart <- function(cluster_genes, top = 3) {
  if (!is.data.frame(cluster_genes) ||
      !all(c("cluster", "gene", "count") %in% colnames(cluster_genes))) {
    stop("`cluster_genes` must be the output of `cluster_top_genes()`.", call. = FALSE)
  }

  cluster_genes <- stats::na.omit(cluster_genes)

  df_summary <- dplyr::group_by(cluster_genes, .data$cluster)
  df_summary <- dplyr::mutate(
    df_summary,
    total_count = sum(.data$count),
    percentage = .data$count / .data$total_count * 100
  )
  df_summary <- dplyr::arrange(df_summary, .data$cluster,
                               dplyr::desc(.data$count))
  df_summary <- dplyr::slice_head(df_summary, n = top)
  df_summary <- dplyr::mutate(
    df_summary,
    gene_label = paste0(.data$gene, " (",
                        round(.data$percentage, 1), "%)")
  )
  df_summary <- as.data.frame(df_summary)

  for (cl in unique(df_summary$cluster)) {
    df_one <- df_summary[df_summary$cluster == cl, , drop = FALSE]
    rest <- max(df_one$total_count) - sum(df_one$count)
    rest_percent <- rest / max(df_one$total_count) * 100

    other_row <- data.frame(
      cluster = cl,
      gene = "other",
      count = rest,
      total_count = max(df_one$total_count),
      percentage = rest_percent,
      gene_label = paste0("other:", round(rest_percent / 100, 2)),
      row.names = NULL
    )
    df_summary <- rbind(df_summary, other_row)
  }

  n_colors <- length(unique(df_summary$gene))
  fill_colors <- grDevices::colorRampPalette(
    ggsci::pal_npg()(10)
  )(n_colors)

  ggplot2::ggplot(df_summary, ggplot2::aes(x = "", y = .data$percentage,
                                           fill = .data$gene)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::facet_wrap(~ paste("Cluster", .data$cluster), nrow = 1) +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(label = .data$gene_label),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.5, color = "white", fontface = "bold") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    ggplot2::labs(title = paste("Top", top, "Genes Percentage by Cluster")) +
    # Stacked percentages can exceed 100 by floating-point roundoff.
    # Clamp the boundary instead of censoring it (which drops a whole slice).
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      oob = function(x, range) pmin(pmax(x, range[1]), range[2])
    ) +
    ggplot2::scale_fill_manual(values = fill_colors)
}
