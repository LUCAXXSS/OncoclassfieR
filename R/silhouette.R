# Silhouette analysis of cluster stability ------------------------------

#' Average silhouette width across k
#'
#' Computes the average silhouette width for a range of `k` (number of
#' clusters) using the same Euclidean distance on the Jaccard similarity
#' representation as [jaccard_cluster()] (approach A: each sample is a point
#' in the space of its similarities to all other samples).  For each `k` the
#' samples are clustered with [stats::hclust()] (ward.D2 by default) and the
#' mean silhouette width is returned, so the trade-off between the number of
#' clusters and the separation of the resulting groups can be inspected.
#'
#' @param cluster_res Either the result list of [jaccard_cluster()] (which
#'   contains a `similarity_matrix` element) or a numeric similarity matrix
#'   (square, samples in rows and columns).
#' @param k_range Integer vector of cluster numbers to evaluate (default
#'   `2:10`). Every value must be `>= 2` and `<=` the number of samples.
#' @param clustering_method Character. Hierarchical clustering method
#'   (default `"ward.D2"`), passed to [stats::hclust()].
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item `silhouette`: data frame with columns `k` and `silhouette`
#'       (mean silhouette width for that `k`).
#'     \item `plot`: a `ggplot` of silhouette width against `k`; call
#'       `print(out$plot)` (or plot it interactively) to draw.
#'   }
#'
#' @export
#' @examples
#' res <- jaccard_cluster(example_maf, n_clusters = 4, plot = FALSE)
#' out <- silhouette_analysis(res, k_range = 2:4)
#' out$silhouette
#' out$plot
silhouette_analysis <- function(cluster_res, k_range = 2:10,
                                clustering_method = "ward.D2") {
  sim <- .get_similarity_matrix(cluster_res)

  dist_mat <- stats::as.dist(stats::dist(sim))

  silhouette_scores <- vapply(k_range, function(k) {
    clusters <- stats::cutree(
      stats::hclust(dist_mat, method = clustering_method), k = k
    )
    sil <- cluster::silhouette(clusters, dist_mat)
    mean(sil[, 3])   # third column is the silhouette width
  }, numeric(1))

  df_sil <- data.frame(k = k_range, silhouette = silhouette_scores)

  p_sil <- ggplot2::ggplot(df_sil, ggplot2::aes(x = .data$k, y = .data$silhouette)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1.2) +
    ggplot2::geom_point(color = "red", size = 3) +
    ggplot2::labs(title = "Silhouette analysis",
                  x = "K",
                  y = "average Silhouette score") +
    ggplot2::theme_classic()

  list(silhouette = df_sil, plot = p_sil)
}

# Extract a similarity matrix from a jaccard_cluster() result list, or
# pass a plain matrix through unchanged.
.get_similarity_matrix <- function(x) {
  if (is.list(x) && !is.data.frame(x) && "similarity_matrix" %in% names(x)) {
    x <- x$similarity_matrix
  }
  if (!is.matrix(x)) {
    stop("`cluster_res` must be a jaccard_cluster() result or a similarity matrix.",
         call. = FALSE)
  }
  x
}
