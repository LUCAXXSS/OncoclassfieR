# Similarity heatmaps drawn in a given subgroup order ---------------------

# Small cohort with four well separated gene sets, so that the clustering is
# stable and the tests can rely on the resulting subgroups.
ordered_maf <- data.frame(
  Tumor_Sample_Barcode = rep(paste0("S", 1:8), each = 2),
  Hugo_Symbol = c("A", "B", "A", "B", "A", "B", "A", "C",
                  "D", "E", "D", "E", "D", "E", "D", "F"),
  stringsAsFactors = FALSE
)

# Two subgroups (max) with samples S1-S4 and S5-S8.
group_labels <- c(S1 = "b", S2 = "b", S3 = "a", S4 = "a",
                  S5 = "c", S6 = "c", S7 = "c", S8 = "c")

quiet_heatmap <- function(...) {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  plot_jaccard_heatmap(...)
}

test_that("jaccard_cluster still returns the original three elements", {
  res <- jaccard_cluster(ordered_maf, n_clusters = 2, plot = FALSE,
                         clustering_method = "ward.D2")

  expect_named(res, c("heatmap", "sample_clusters", "similarity_matrix"))
  expect_s3_class(res$heatmap, "pheatmap")
  expect_equal(sort(unique(res$sample_clusters$cluster)), 1:2)
  expect_equal(nrow(res$sample_clusters), 8L)

  # The clustering result is not affected by the plotting helper.
  hm <- quiet_heatmap(res, group_order = c(2, 1))
  expect_equal(res$sample_clusters$Tumor_Sample_Barcode,
               names(stats::cutree(res$heatmap$tree_row, k = 2)))
})

test_that("plot_jaccard_heatmap orders the samples by subgroup", {
  res <- jaccard_cluster(ordered_maf, n_clusters = 2, plot = FALSE,
                         clustering_method = "ward.D2")
  hm <- quiet_heatmap(res, group_order = c("2", "1"))

  expect_s3_class(hm$heatmap, "pheatmap")
  expect_equal(hm$group_order, c("2", "1"))
  expect_equal(as.character(hm$group$cluster), c("2", "2", "2", "2",
                                                 "1", "1", "1", "1"))
  expect_equal(hm$sample_order, hm$group$Tumor_Sample_Barcode)

  # Blocks are contiguous and separated by gaps on both axes.
  expect_equal(unique(as.character(hm$group$cluster)),
               c("2", "1"))
  expect_equal(hm$gaps, 4)
  expect_equal(rownames(hm$annotation), hm$sample_order)
})

test_that("plot_jaccard_heatmap uses the clustering of the input result", {
  res <- jaccard_cluster(ordered_maf, n_clusters = 2, plot = FALSE,
                         clustering_method = "ward.D2")

  # Without group_order the subgroups keep the numeric-aware default order.
  hm_default <- quiet_heatmap(res)
  expect_equal(hm_default$group_order, c("1", "2"))

  # Reversing the order reverses the blocks.
  hm_rev <- quiet_heatmap(res, group_order = c(2, 1))
  expect_setequal(unlist(split(hm_rev$sample_order,
                               as.character(hm_rev$group$cluster))),
                  hm_rev$sample_order)
  expect_equal(unique(as.character(hm_rev$group$cluster)), c("2", "1"))
  expect_setequal(hm_rev$sample_order, hm_default$sample_order)
})

test_that("plot_jaccard_heatmap accepts a plain similarity matrix", {
  sim <- jaccard_similarity_matrix(ordered_maf)
  hm <- quiet_heatmap(sim, group = group_labels, group_order = c("c", "a", "b"))

  expect_equal(hm$group_order, c("c", "a", "b"))
  expect_equal(as.character(hm$group$cluster),
               c(rep("c", 4), rep("a", 2), rep("b", 2)))
  expect_setequal(hm$sample_order[1:4], c("S5", "S6", "S7", "S8"))
  expect_setequal(hm$sample_order[5:6], c("S3", "S4"))
  expect_setequal(hm$sample_order[7:8], c("S1", "S2"))
  expect_equal(hm$gaps, c(4, 6))

  # A matrix without a grouping cannot be drawn.
  expect_error(quiet_heatmap(sim), "group")
})

test_that("ordering modes and grouping shapes agree", {
  sim <- jaccard_similarity_matrix(ordered_maf)

  by_input <- quiet_heatmap(sim, group = group_labels,
                            group_order = c("a", "b", "c"),
                            order_within_group = "input")
  by_sim <- quiet_heatmap(sim, group = group_labels,
                          order_within_group = "similarity")

  # "input" keeps the order of the grouping vector inside each block.
  expect_equal(by_input$sample_order,
               c("S3", "S4", "S1", "S2", "S5", "S6", "S7", "S8"))
  # Default group order is the numeric-aware sort of the labels.
  expect_equal(by_sim$group_order, c("a", "b", "c"))

  # A named list and a data frame describe the same grouping.
  as_list <- list(a = c("S3", "S4"), b = c("S1", "S2"),
                  c = c("S5", "S6", "S7", "S8"))
  as_df <- data.frame(Tumor_Sample_Barcode = names(group_labels),
                      subgroup = unname(group_labels),
                      stringsAsFactors = FALSE)

  hm_list <- quiet_heatmap(sim, group = as_list)
  hm_df <- quiet_heatmap(sim, group = as_df, group_order = c("a", "b", "c"))

  expect_equal(hm_list$sample_order, hm_df$sample_order)
  expect_equal(names(hm_df$annotation), "subgroup")
})

test_that("plot_jaccard_heatmap handles colours and extra annotations", {
  sim <- jaccard_similarity_matrix(ordered_maf)
  ann <- data.frame(cohort = rep(c("p", "q"), 4), row.names = rownames(sim),
                    stringsAsFactors = FALSE)

  hm <- quiet_heatmap(sim, group = group_labels,
                      group_order = c("c", "b", "a"),
                      group_colors = c(c = "red"), annotation = ann)

  expect_equal(colnames(hm$annotation), c("Subgroup", "cohort"))
  expect_equal(levels(hm$annotation$Subgroup), c("c", "b", "a"))
  expect_equal(as.character(hm$annotation$Subgroup[1]), "c")
  expect_equal(rownames(hm$annotation), hm$sample_order)
})

test_that("samples missing from the grouping are drawn last with a warning", {
  sim <- jaccard_similarity_matrix(ordered_maf)

  expect_warning(
    hm <- quiet_heatmap(sim, group = c(S1 = "a", S2 = "a")),
    "missing from `group`"
  )

  expect_setequal(hm$group_order, c("a", "Not assigned"))
  expect_setequal(hm$sample_order[1:2], c("S1", "S2"))
  expect_setequal(tail(hm$sample_order, 6),
                  c("S3", "S4", "S5", "S6", "S7", "S8"))
})

test_that("unknown subgroup labels in group_order are reported", {
  sim <- jaccard_similarity_matrix(ordered_maf)

  expect_warning(
    hm <- quiet_heatmap(sim, group = group_labels, group_order = c("z", "b")),
    "absent from `group`"
  )

  expect_equal(hm$group_order, c("b", "a", "c"))
})

test_that("a factor grouping keeps its level order", {
  sim <- jaccard_similarity_matrix(ordered_maf)
  grp <- factor(c("late", "late", "early", "early", "mid", "mid", "mid", "mid"),
                levels = c("early", "mid", "late"))
  names(grp) <- paste0("S", 1:8)

  hm <- quiet_heatmap(sim, group = grp)

  expect_equal(hm$group_order, c("early", "mid", "late"))
  expect_setequal(tail(hm$sample_order, 2), c("S1", "S2"))
})

test_that("plot_jaccard_heatmap can return the plot without drawing it", {
  res <- jaccard_cluster(ordered_maf, n_clusters = 2, plot = FALSE,
                         clustering_method = "ward.D2")

  # plot = FALSE builds the pheatmap object silently.
  hm <- plot_jaccard_heatmap(res, group_order = c(2, 1), plot = FALSE)

  expect_s3_class(hm$heatmap, "pheatmap")
  expect_equal(hm$group_order, c("2", "1"))
  expect_equal(hm$sample_order, hm$group$Tumor_Sample_Barcode)
})
