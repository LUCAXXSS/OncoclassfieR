# Core Jaccard similarity clustering -------------------------------------

#' Keep the most frequently mutated genes in a MAF object
#'
#' Reduces a MAF object to the `max_features` most frequently mutated genes.
#' This is the pre-processing step used before Jaccard clustering so that
#' rare genes do not dominate the pairwise sample similarities.
#'
#' @param maf An object of class `MAF` (from the `maftools` package), or any
#'   object with a `@data` slot containing a MAF-like data frame with a
#'   `Hugo_Symbol` column. A plain data frame with a `Hugo_Symbol` column is
#'   also accepted and returned unchanged apart from the gene filtering.
#' @param max_features Integer. Number of top genes to keep.
#' @param drop_empty_samples Logical. When `FALSE` (the default) samples that
#'   have no mutation among the retained genes are kept in the analysis: the
#'   returned data frame carries the full sample universe in
#'   `attr(out, "samples")`, so [jaccard_similarity_matrix()] and
#'   [jaccard_cluster()] still score them (as all-zero similarity profiles).
#'   This reproduces the behaviour of the original analysis pipeline. When
#'   `TRUE`, only samples with at least one mutation among the retained genes
#'   enter the similarity matrix.
#'
#' @return A data frame (the content of `maf@data`) restricted to the
#'   `max_features` most frequent genes. The attribute `samples` holds the
#'   samples that should enter the similarity matrix (see
#'   `drop_empty_samples`); the rows of the data frame only cover the samples
#'   that actually carry one of the retained genes.
#'
#' @export
#' @examples
#' \donttest{
#' maf_top <- preprocess_maf_top_genes(example_maf, max_features = 20)
#' length(attr(maf_top, "samples"))
#'
#' # Strict version: drop samples without any mutation in the top genes
#' maf_top_strict <- preprocess_maf_top_genes(example_maf, max_features = 20,
#'                                            drop_empty_samples = TRUE)
#' }
preprocess_maf_top_genes <- function(maf, max_features = 50,
                                     drop_empty_samples = FALSE) {
  # The sample universe has to be read before the gene filter reduces the
  # rows; a factor `Tumor_Sample_Barcode` defines it in level order, exactly
  # like `names(table(...))` did in the original pipeline.
  samples <- .sample_universe(maf)

  maf_data <- .maf_to_dataframe(maf)

  gene_counts <- table(maf_data$Hugo_Symbol)
  gene_counts <- sort(gene_counts, decreasing = TRUE)

  if (length(gene_counts) == 0) {
    stop("No genes found in the provided MAF data.", call. = FALSE)
  }

  top_genes <- names(gene_counts)[seq_len(min(max_features, length(gene_counts)))]
  out <- maf_data[maf_data$Hugo_Symbol %in% top_genes, , drop = FALSE]

  if (isTRUE(drop_empty_samples)) {
    samples <- unique(as.character(out$Tumor_Sample_Barcode))
  }

  attr(out, "samples") <- samples
  out
}

#' Compute the pairwise Jaccard similarity matrix of samples
#'
#' Given a MAF-like data frame with `Tumor_Sample_Barcode` and `Hugo_Symbol`
#' columns, computes the Jaccard similarity (size of gene-set intersection /
#' size of gene-set union) between every pair of samples. Duplicated
#' sample-gene pairs are removed first.
#'
#' The diagonal (self-similarity, always 1) is replaced by the largest
#' observed off-diagonal similarity, and the whole matrix is raised to the
#' power `enhancement`, mirroring the original pipeline so that the heatmap
#' colours are compressed towards the low end.
#'
#' Samples without any mutation among the genes present in `df` score
#' `0 / 0 = NaN` against each other; these entries are set to 0, so such
#' samples appear as all-zero rows (and columns). This reproduces the original
#' pipeline, where the sample universe was taken from the factor levels of
#' `Tumor_Sample_Barcode` and empty samples were kept as zero profiles. Use
#' `drop_empty_samples = TRUE` in [preprocess_maf_top_genes()] to exclude them
#' instead.
#'
#' @param df Data frame (or `MAF` object) containing at least
#'   `Tumor_Sample_Barcode` and `Hugo_Symbol` columns. The sample universe is
#'   taken from `attr(df, "samples")` when present (as returned by
#'   [preprocess_maf_top_genes()]) and otherwise from the samples observed in
#'   `df`.
#' @param enhancement Numeric between 0 and 1. Exponent applied to the
#'   similarity matrix before plotting/clustering (default 0.6).
#'
#' @return A numeric matrix with samples in both row and column names. Samples
#'   that carry no mutation among the genes of `df` are included as all-zero
#'   rows unless they were dropped upstream.
#'
#' @export
#' @examples
#' sim <- jaccard_similarity_matrix(example_maf)
#' dim(sim)
jaccard_similarity_matrix <- function(df, enhancement = 0.6) {
  # Read the sample universe before the de-duplication step, which only keeps
  # the `Tumor_Sample_Barcode` / `Hugo_Symbol` columns.
  samples <- .similarity_samples(df)

  df <- .maf_to_dataframe(df)
  df <- .dedup_sample_gene(df)

  if (length(samples) == 0) {
    stop("No samples found in the provided MAF data.", call. = FALSE)
  }

  n <- length(samples)

  sim <- matrix(0, nrow = n, ncol = n,
                dimnames = list(samples, samples))

  # Gene sets of the samples that survived the gene filter; samples kept only
  # through `attr(df, "samples")` are absent and handled as empty sets below.
  gene_sets <- split(df$Hugo_Symbol, df$Tumor_Sample_Barcode)

  for (i in samples) {
    genes_i <- gene_sets[[i]]
    if (is.null(genes_i)) genes_i <- character(0)
    for (j in samples) {
      genes_j <- gene_sets[[j]]
      if (is.null(genes_j)) genes_j <- character(0)
      sim[i, j] <- .jaccard(genes_i, genes_j)
    }
  }

  # Empty gene sets give 0 / 0 = NaN against each other; the original pipeline
  # scored those pairs 0. This has to happen before `max(sim)` below: a NaN
  # would otherwise propagate into the self-similarity replacement and make
  # the whole matrix NaN, breaking the downstream clustering.
  sim[is.nan(sim)] <- 0

  # Self-similarity (1) would dominate the heatmap / clustering; replace it
  # with the largest observed off-diagonal similarity as in the pipeline.
  sim[sim == 1] <- -1
  sim[sim == -1] <- max(sim)

  sim <- sim^enhancement
  sim
}

#' Cluster samples by Jaccard similarity of their mutation profiles
#'
#' Runs the full Jaccard clustering workflow: pre-processes the MAF data
#' (optional), computes the pairwise Jaccard similarity matrix, draws a
#' pheatmap with row/column hierarchical clustering, and cuts the row tree
#' into `n_clusters` sample groups.
#'
#' @param df Data frame (or `MAF` object) with `Tumor_Sample_Barcode` and
#'   `Hugo_Symbol` columns.
#' @param n_clusters Integer. Number of sample clusters (passed to
#'   [stats::cutree()]).
#' @param clustering_method Character. Hierarchical clustering method used by
#'   [pheatmap::pheatmap()] (e.g. `"complete"`, `"ward.D2"`).
#' @param enhancement Numeric between 0 and 1. Exponent applied to the
#'   similarity matrix (see [jaccard_similarity_matrix()]).
#' @param show_colnames Logical. Whether to show sample names in the heatmap.
#' @param plot Logical. Whether to draw the heatmap on the current graphics
#'   device.
#' @param ... Further arguments passed to [pheatmap::pheatmap()].
#'
#' @return A list with three elements:
#'   \itemize{
#'     \item `heatmap`: the pheatmap object (contains `tree_row` used for
#'       the clustering).
#'     \item `sample_clusters`: data frame with columns
#'       `Tumor_Sample_Barcode` and `cluster`.
#'     \item `similarity_matrix`: the enhanced Jaccard similarity matrix.
#'   }
#'
#' @details
#' The samples that enter the clustering are the sample universe of `df` (see
#' [jaccard_similarity_matrix()]); when `df` comes from
#' [preprocess_maf_top_genes()] this includes samples without any mutation
#' among the retained genes, which are clustered as all-zero similarity
#' profiles. Use `drop_empty_samples = TRUE` there to exclude them.
#'
#' The heatmap is drawn in the order of the hierarchical clustering. To draw
#' the same similarity matrix with a known grouping, in a hand-picked subgroup
#' order, use [plot_jaccard_heatmap()].
#'
#' @seealso [plot_jaccard_heatmap()] to re-draw the similarity matrix with a
#'   given subgroup order, and [silhouette_analysis()] to choose `n_clusters`.
#'
#' @export
#' @examples
#' \donttest{
#' res <- jaccard_cluster(example_maf, n_clusters = 3, plot = FALSE)
#' head(res$sample_clusters)
#' }
jaccard_cluster <- function(df,
                            n_clusters,
                            clustering_method = "complete",
                            enhancement = 0.6,
                            show_colnames = TRUE,
                            plot = TRUE,
                            ...) {
  df <- .maf_to_dataframe(df)
  df <- .dedup_sample_gene(df)

  sim <- jaccard_similarity_matrix(df, enhancement = enhancement)

  if (plot) {
    p <- pheatmap::pheatmap(
      sim,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      scale = "none",
      show_rownames = FALSE,
      show_colnames = show_colnames,
      color = viridis::magma(100)[20:100],
      clustering_method = clustering_method,
      cutree_rows = n_clusters,
      cutree_cols = n_clusters,
      border_color = NA,
      ...
    )
  } else {
    # Still build the clustering trees without drawing anything.
    p <- pheatmap::pheatmap(
      sim,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      scale = "none",
      show_rownames = FALSE,
      show_colnames = show_colnames,
      color = viridis::magma(100)[20:100],
      clustering_method = clustering_method,
      cutree_rows = n_clusters,
      cutree_cols = n_clusters,
      border_color = NA,
      silent = TRUE,
      ...
    )
  }

  sample_clusters <- stats::cutree(p$tree_row, k = n_clusters)
  sample_clusters <- data.frame(
    Tumor_Sample_Barcode = names(sample_clusters),
    cluster = unname(sample_clusters),
    row.names = NULL
  )

  list(
    heatmap = p,
    sample_clusters = sample_clusters,
    similarity_matrix = sim
  )
}

# Internal helpers --------------------------------------------------------

.maf_to_dataframe <- function(x) {
  sample_universe <- attr(x, "samples", exact = TRUE)

  if (methods::is(x, "MAF")) {
    x <- x@data
  } else if (is.data.frame(x)) {
    x <- x
  } else {
    stop("`x` must be a MAF object or a data frame.", call. = FALSE)
  }
  # MAF objects store their data as a data.table, which breaks the base
  # `df[, col]` column selection used below; convert to a plain data frame.
  # Factor columns are also converted to character so that unused factor
  # levels (e.g. the full sample set of a merged MAF) do not leak into
  # `table()`, `%in%` and `split()` downstream.
  if (methods::is(x, "data.table")) {
    x <- as.data.frame(x)
  }
  if (is.data.frame(x)) {
    factor_cols <- vapply(x, is.factor, logical(1))
    x[factor_cols] <- lapply(x[factor_cols], as.character)
  }

  # Preserve the sample universe carried by preprocess_maf_top_genes().
  if (!is.null(sample_universe)) {
    attr(x, "samples") <- sample_universe
  }
  x
}

# Sample universe of a MAF object / data frame: the level order of a factor
# `Tumor_Sample_Barcode` (the semantics of `names(table(...))` in the original
# pipeline), restricted to the samples that are actually present, otherwise
# the samples in order of first appearance.
.sample_universe <- function(x) {
  if (methods::is(x, "MAF")) {
    barcode <- x@data[["Tumor_Sample_Barcode"]]
  } else if (is.data.frame(x)) {
    barcode <- x[["Tumor_Sample_Barcode"]]
  } else {
    stop("`x` must be a MAF object or a data frame.", call. = FALSE)
  }

  if (is.null(barcode)) {
    return(character(0))
  }

  observed <- unique(as.character(barcode))

  if (is.factor(barcode)) {
    levels_present <- levels(barcode)
    return(levels_present[levels_present %in% observed])
  }

  observed
}

# Samples to put into the similarity matrix: the universe attached by
# preprocess_maf_top_genes() when available, otherwise the observed samples.
.similarity_samples <- function(df) {
  samples <- attr(df, "samples", exact = TRUE)

  if (!is.null(samples)) {
    return(unique(as.character(samples)))
  }

  .sample_universe(df)
}

.dedup_sample_gene <- function(df) {
  needed <- c("Tumor_Sample_Barcode", "Hugo_Symbol")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop("Input data must contain column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  sample_universe <- attr(df, "samples", exact = TRUE)

  df <- df[, needed]
  df$check <- paste0(df$Tumor_Sample_Barcode, df$Hugo_Symbol)
  out <- df[!duplicated(df$check), needed, drop = FALSE]

  if (!is.null(sample_universe)) {
    attr(out, "samples") <- sample_universe
  }
  out
}

.jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}
