# Similarity heatmaps in a given subgroup order ----------------------------
#
# jaccard_cluster() draws the heatmap in the order of the hierarchical
# clustering. A figure often has to follow an existing subgroup assignment
# instead, with the subgroup blocks in a hand-picked order (so that the panel
# lines up with the oncoplots, survival curves, ...). plot_jaccard_heatmap()
# does exactly that: it takes the similarity matrix (or the result of
# jaccard_cluster()) and re-draws it with the samples arranged in subgroup
# blocks, the blocks following `group_order`.

#' Draw a Jaccard similarity heatmap in a given subgroup order
#'
#' Re-draws the pairwise Jaccard similarity matrix as a `pheatmap` in which the
#' samples are arranged in subgroup blocks instead of by hierarchical
#' clustering. The subgroups themselves follow `group_order`, so a figure can
#' be reproduced with the same block order as the other panels.
#'
#' @param x The result list of [jaccard_cluster()], or a similarity matrix as
#'   returned by [jaccard_similarity_matrix()] (square, unique sample names in
#'   the row names).
#' @param group Sample grouping, in any of these forms:
#'   \itemize{
#'     \item a data frame with a sample column (`Tumor_Sample_Barcode`,
#'       `sample`, `Sample`, `ID`, `Barcode`) and a subgroup column (`cluster`,
#'       `group`, `subgroup`, `subtype`); with two unnamed columns the first is
#'       taken as the sample and the second as the subgroup. A factor subgroup
#'       column keeps its level order as the default display order,
#'     \item a named vector, `c(S1 = "1", S2 = "3", S3 = "1")`, whose names are
#'       the sample barcodes and whose values are the subgroup labels,
#'     \item a named list of sample vectors, one element per subgroup,
#'     \item `NULL` (the default): the `sample_clusters` element of `x` is used
#'       when `x` is a [jaccard_cluster()] result. For a plain matrix `group`
#'       is required.
#'   }
#' @param group_order Character (or numeric) vector with the subgroups in the
#'   order in which their blocks should appear in the heatmap and in the
#'   annotation legend, e.g. `c(3, 1, 4, 2)`. Subgroups that are present in
#'   `group` but absent from `group_order` are appended at the end, in order of
#'   first appearance; labels that are absent from `group` are ignored (with a
#'   warning). When `NULL`, the factor levels of a factor grouping are used,
#'   otherwise a numeric-aware sort of the subgroup labels.
#' @param group_colors Named character vector of colours, one per subgroup
#'   (`c("1" = "#E64B35", "2" = "#4DBBD5", ...)`). Subgroups without a colour
#'   receive one from the default palette (NPG colours).
#' @param annotation Optional data frame of extra sample annotations, with the
#'   sample barcodes as `rownames`, drawn as additional annotation tracks next
#'   to the subgroup bar. Use `annotation_colors` (passed on to
#'   [pheatmap::pheatmap()]) to set their colours.
#' @param order_within_group How the samples inside one subgroup are ordered:
#'   `"hclust"` (default) hierarchical clustering of the similarity submatrix
#'   with `clustering_method`, `"similarity"` by decreasing mean similarity to
#'   the other samples of the subgroup, or `"input"` by their order in `group`
#'   (for samples that are not in `group`, the order of the matrix).
#' @param clustering_method Character. Clustering method used for the
#'   within-subgroup ordering (`"complete"`, `"ward.D2"`, ...).
#' @param show_colnames Logical. Whether to show sample names on the columns.
#' @param show_rownames Logical. Whether to show sample names on the rows.
#' @param color Heatmap colour vector (default: the `viridis::magma` ramp used
#'   by [jaccard_cluster()]).
#' @param plot Logical. Whether to draw the heatmap on the current graphics
#'   device. `FALSE` still builds and returns the `pheatmap` object.
#' @param ... Further arguments passed to [pheatmap::pheatmap()] (for example
#'   `main`, `fontsize`, `annotation_colors`, `legend`).
#'
#' @return A list with the elements
#'   \itemize{
#'     \item `heatmap`: the `pheatmap` object; print it (or use
#'       `$gtable`) to draw the figure again.
#'     \item `sample_order`: character vector of the samples in the order in
#'       which they are drawn (matrix rows and columns alike).
#'     \item `group`: data frame with columns `Tumor_Sample_Barcode` and
#'       `cluster` (the subgroup labels), sorted like the heatmap.
#'     \item `group_order`: the subgroup order actually used (after appending
#'       the subgroups missing from `group_order`).
#'     \item `annotation`: the annotation data frame drawn next to the heatmap.
#'     \item `gaps`: the row/column positions at which the block gaps are drawn.
#'   }
#'
#' @details
#' The heatmap of [jaccard_cluster()] follows the clustering order. This
#' function leaves that result untouched and draws the same similarity matrix
#' again, ordered by the subgroups:
#'
#' ```r
#' res <- jaccard_cluster(maf_top, n_clusters = 4, plot = FALSE)
#' hm <- plot_jaccard_heatmap(res, group_order = c(3, 1, 4, 2))
#' hm$sample_order     # the samples in the order they are drawn
#' ```
#'
#' Samples of the similarity matrix that are missing from `group` are drawn
#' last as a `"Not assigned"` block (with a warning) rather than dropped.
#' Samples of `group` that are not part of the matrix are ignored.
#'
#' @seealso [jaccard_cluster()] for the clustering itself,
#'   [cluster_top_genes()], [plot_cluster_oncoplot()] and
#'   [plot_tmb_by_cluster()] for per-subgroup analyses that can consume
#'   `hm$group`.
#'
#' @export
#' @examples
#' \donttest{
#' res <- jaccard_cluster(example_maf, n_clusters = 4, plot = FALSE)
#' hm <- plot_jaccard_heatmap(res, group_order = c(4, 2, 1, 3))
#' hm$sample_order
#'
#' # Custom colours, another order and a table of extra annotations
#' cohort <- data.frame(cohort = rep(c("A", "B"), length.out = ncol(res$similarity_matrix)),
#'                      row.names = colnames(res$similarity_matrix))
#' hm <- plot_jaccard_heatmap(res, group_order = c(2, 4), 
#'                            group_colors = c("2" = "#E64B35"),
#'                            annotation = cohort)
#'
#' # A hand-made grouping also works
#' my_group <- data.frame(
#'   Tumor_Sample_Barcode = res$sample_clusters$Tumor_Sample_Barcode,
#'   subgroup = res$sample_clusters$cluster
#' )
#' plot_jaccard_heatmap(res$similarity_matrix, group = my_group,
#'                      group_order = c(2, 4), order_within_group = "input")
#' }
plot_jaccard_heatmap <- function(x,
                                 group = NULL,
                                 group_order = NULL,
                                 group_colors = NULL,
                                 annotation = NULL,
                                 order_within_group = c("hclust", "similarity",
                                                        "input"),
                                 clustering_method = "complete",
                                 show_colnames = TRUE,
                                 show_rownames = FALSE,
                                 color = NULL,
                                 plot = TRUE,
                                 ...) {
  order_within_group <- match.arg(order_within_group)

  is_result <- is.list(x) && !is.data.frame(x) &&
    "similarity_matrix" %in% names(x)
  if (is_result && is.null(group) && "sample_clusters" %in% names(x)) {
    group <- x$sample_clusters
  }
  if (is.null(group)) {
    stop("`group` is required when `x` is a similarity matrix; pass a data ",
         "frame, named vector or named list of subgroup labels. It can only ",
         "default to the grouping when `x` is a `jaccard_cluster()` result.",
         call. = FALSE)
  }

  sim <- .get_similarity_matrix(x)

  samples <- rownames(sim)
  if (is.null(samples) || anyDuplicated(samples)) {
    stop("`x` must be a similarity matrix with unique sample row names.",
         call. = FALSE)
  }

  dots <- list(...)
  conflicting <- intersect(names(dots),
                           c("annotation", "annotation_col", "annotation_row"))
  if (length(conflicting) > 0) {
    stop("Extra sample annotations must be passed through the `annotation` ",
         "argument, not through `...` (conflicting argument(s): ",
         paste(conflicting, collapse = ", "), ").", call. = FALSE)
  }

  if (is.null(color)) {
    color <- viridis::magma(100)[20:100]
  }

  grouping <- .normalise_group(group, samples, group_order)

  .plot_grouped_heatmap(
    sim, grouping = grouping, group_order = group_order,
    group_colors = group_colors, annotation = annotation,
    order_within_group = order_within_group,
    clustering_method = clustering_method,
    show_colnames = show_colnames, show_rownames = show_rownames,
    color = color, plot = plot, dots = dots
  )
}

# Internal: heatmap ordered by the subgroup blocks.
.plot_grouped_heatmap <- function(sim, grouping, group_order = NULL,
                                  group_colors = NULL, annotation = NULL,
                                  order_within_group = "hclust",
                                  clustering_method = "complete",
                                  show_colnames = TRUE, show_rownames = FALSE,
                                  color, plot = TRUE, dots = list()) {
  samples <- rownames(sim)
  group_name <- attr(grouping, "group_name")

  # ---- subgroup order --------------------------------------------------
  present <- unique(grouping$group)
  if (is.null(group_order)) {
    group_order <- attr(grouping, "group_order")
  }
  if (is.null(group_order)) {
    group_order <- .natural_sort(present)
  }
  group_order <- as.character(group_order)

  appended <- setdiff(present, group_order)
  if (length(appended) > 0) {
    group_order <- c(group_order, appended)
  }
  unused <- setdiff(group_order, present)
  if (length(unused) > 0) {
    warning("`group_order` contains subgroup(s) absent from `group`: ",
            paste(unused, collapse = ", "),
            ". Check the subgroup labels.", call. = FALSE)
    group_order <- setdiff(group_order, unused)
  }

  # ---- align the grouping with the matrix ------------------------------
  idx <- match(samples, grouping$sample)
  if (anyNA(idx)) {
    unassigned <- samples[is.na(idx)]
    warning(length(unassigned), " sample(s) are missing from `group` and are ",
            "drawn last as subgroup \"Not assigned\": ",
            paste(utils::head(unassigned, 5), collapse = ", "),
            if (length(unassigned) > 5) ", ..." else "", call. = FALSE)
    grouping <- rbind(
      grouping,
      data.frame(sample = unassigned, group = "Not assigned",
                 order_in = seq_along(unassigned), stringsAsFactors = FALSE)
    )
    idx <- match(samples, grouping$sample)
    group_order <- c(group_order, setdiff("Not assigned", group_order))
  }

  sample_group <- grouping$group[idx]
  sample_rank <- grouping$order_in[idx]

  # ---- order the samples: subgroup blocks, order inside the block ------
  ord <- integer(0)
  for (g in group_order) {
    pos <- which(sample_group == g)
    if (length(pos) == 0) next

    if (length(pos) > 1 && identical(order_within_group, "hclust")) {
      sub <- sim[pos, pos, drop = FALSE]
      pos <- pos[stats::hclust(stats::dist(sub),
                               method = clustering_method)$order]
    } else if (length(pos) > 1 && identical(order_within_group, "similarity")) {
      sub <- sim[pos, pos, drop = FALSE]
      score <- (rowSums(sub) - diag(sub)) / (length(pos) - 1)
      pos <- pos[order(-score)]
    } else if (identical(order_within_group, "input")) {
      pos <- pos[order(sample_rank[pos])]
    }

    ord <- c(ord, pos)
  }

  mat <- sim[ord, ord, drop = FALSE]
  blocks <- factor(sample_group[ord], levels = group_order)

  # ---- annotation table ------------------------------------------------
  ann <- data.frame(blocks, row.names = samples[ord], stringsAsFactors = FALSE)
  names(ann) <- group_name

  if (!is.null(annotation)) {
    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
    if (is.null(rownames(annotation))) {
      stop("`annotation` needs the sample barcodes as rownames.", call. = FALSE)
    }
    missing_rows <- setdiff(samples[ord], rownames(annotation))
    if (length(missing_rows) > 0) {
      stop("`annotation` has no row for ", length(missing_rows),
           " sample(s) of the heatmap, e.g. ",
           paste(utils::head(missing_rows, 3), collapse = ", "), ".",
           call. = FALSE)
    }
    ann <- cbind(ann, annotation[samples[ord], , drop = FALSE])
  }

  # ---- colours ---------------------------------------------------------
  group_colors <- .resolve_group_colors(group_order, group_colors)
  ann_colors <- stats::setNames(list(group_colors), group_name)
  if (!is.null(dots$annotation_colors)) {
    ann_colors <- c(ann_colors, dots$annotation_colors)
    dots$annotation_colors <- NULL
  }

  # ---- gaps between the subgroup blocks --------------------------------
  counts <- as.integer(table(blocks))
  counts <- counts[counts > 0]
  gaps <- if (length(counts) > 1) utils::head(cumsum(counts), -1) else NULL

  args <- list(
    mat = mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = color,
    annotation_row = ann,
    annotation_col = ann,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    border_color = NA,
    gaps_row = gaps,
    gaps_col = gaps,
    silent = !isTRUE(plot)
  )
  args[names(dots)] <- dots

  p <- do.call(pheatmap::pheatmap, args)

  list(
    heatmap = p,
    sample_order = samples[ord],
    group = data.frame(Tumor_Sample_Barcode = samples[ord],
                       cluster = as.character(blocks),
                       stringsAsFactors = FALSE),
    group_order = group_order,
    annotation = ann,
    gaps = gaps
  )
}

# Internal: normalise the many accepted `group` shapes into a two-column
# data frame (`sample`, `group`) carrying the subgroup name and the default
# display order as attributes.
.normalise_group <- function(group, samples, group_order = NULL) {
  if (is.null(group)) {
    return(NULL)
  }

  group_name <- "Subgroup"

  if (is.data.frame(group)) {
    if (ncol(group) < 2) {
      stop("`group` needs at least two columns: the sample barcode and the ",
           "subgroup label.", call. = FALSE)
    }
    sample_col <- .match_column(
      group, c("Tumor_Sample_Barcode", "sample", "Sample", "SAMPLES",
               "samples", "ID", "Barcode", "barcode")
    )
    group_col <- .match_column(
      group, c("cluster", "Cluster", "group", "Group", "subgroup",
               "Subgroup", "subtype", "Subtype")
    )
    if (is.null(sample_col) || is.null(group_col) ||
        identical(sample_col, group_col)) {
      sample_col <- colnames(group)[1]
      group_col <- colnames(group)[2]
    }

    values <- group[[group_col]]
    if (is.null(group_order) && is.factor(values)) {
      group_order <- levels(droplevels(values))
    }
    out <- data.frame(
      sample = as.character(group[[sample_col]]),
      group = as.character(values),
      stringsAsFactors = FALSE
    )
    group_name <- group_col
  } else if (is.list(group)) {
    if (is.null(names(group)) || !all(nzchar(names(group)))) {
      stop("When `group` is a list it must be named, one element per subgroup.",
           call. = FALSE)
    }
    out <- do.call(rbind, lapply(names(group), function(g) {
      data.frame(sample = as.character(group[[g]]), group = g,
                 stringsAsFactors = FALSE)
    }))
    if (is.null(group_order)) {
      group_order <- names(group)
    }
  } else if (is.atomic(group)) {
    if (is.null(names(group)) || !all(nzchar(names(group)))) {
      stop("When `group` is a vector it must be named, with the sample ",
           "barcodes as names.", call. = FALSE)
    }
    if (is.null(group_order) && is.factor(group)) {
      group_order <- levels(droplevels(group))
    }
    out <- data.frame(sample = names(group), group = as.character(group),
                      stringsAsFactors = FALSE)
  } else {
    stop("`group` must be a data frame, a named vector, a named list or NULL.",
         call. = FALSE)
  }

  if (nrow(out) == 0) {
    stop("`group` does not contain any sample.", call. = FALSE)
  }

  # A sample listed twice keeps its first assignment; the row order of `group`
  # is also the "input" order used inside a subgroup.
  dup <- duplicated(out$sample)
  if (any(dup)) {
    warning(sum(dup), " sample(s) are listed more than once in `group`; ",
            "keeping the first assignment.", call. = FALSE)
    out <- out[!dup, , drop = FALSE]
  }

  out$order_in <- seq_len(nrow(out))
  attr(out, "group_name") <- group_name
  attr(out, "group_order") <- group_order
  out
}

# Internal: first matching column name.
.match_column <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NULL else hit[1]
}

# Internal: numeric-aware sort of the subgroup labels ("1" < "2" < "10").
.natural_sort <- function(x) {
  x <- as.character(x)
  num <- suppressWarnings(as.numeric(x))
  if (length(num) > 0 && !anyNA(num)) {
    x[order(num)]
  } else {
    sort(x)
  }
}

# Internal: default subgroup colours, completing a partially specified vector.
.resolve_group_colors <- function(group_order, group_colors = NULL) {
  n <- length(group_order)

  if (is.null(group_colors)) {
    out <- .default_group_colors(n)
    names(out) <- group_order
  } else {
    if (is.null(names(group_colors)) || !all(nzchar(names(group_colors)))) {
      stop("`group_colors` must be a named vector, one colour per subgroup.",
           call. = FALSE)
    }
    missing <- setdiff(group_order, names(group_colors))
    if (length(missing) > 0) {
      fill <- .default_group_colors(length(missing))
      names(fill) <- missing
      group_colors <- c(group_colors, fill)
    }
    out <- group_colors[group_order]
  }

  if ("Not assigned" %in% names(out)) {
    out["Not assigned"] <- "#BEBEBE"
  }
  out
}

# Internal: NPG palette, extended by interpolation beyond 10 subgroups.
.default_group_colors <- function(n) {
  if (n <= 0) {
    return(character(0))
  }
  if (n <= 10) {
    ggsci::pal_npg()(n)
  } else {
    grDevices::colorRampPalette(ggsci::pal_npg()(10))(n)
  }
}
