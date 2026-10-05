test_that("preprocess_maf_top_genes keeps only top genes", {
  dat <- panel_subset()
  out <- preprocess_maf_top_genes(dat, max_features = 20)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Tumor_Sample_Barcode", "Hugo_Symbol") %in% colnames(out)))
  expect_lte(length(unique(out$Hugo_Symbol)), 20)
  expect_equal(nrow(out), nrow(dat[dat$Hugo_Symbol %in%
                                     unique(out$Hugo_Symbol), ]))
})

test_that("jaccard_similarity_matrix returns a square matrix", {
  dat <- panel_subset()
  sim <- jaccard_similarity_matrix(dat)
  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), ncol(sim))
  expect_equal(nrow(sim), length(unique(dat$Tumor_Sample_Barcode)))
  expect_equal(rownames(sim), colnames(sim))
  expect_true(all(sim >= 0 & sim <= 1))
})

test_that("jaccard_cluster returns clusters, heatmap and matrix", {
  dat <- panel_subset()
  res <- jaccard_cluster(dat, n_clusters = 4, plot = FALSE)
  expect_type(res, "list")
  expect_named(res, c("heatmap", "sample_clusters", "similarity_matrix"))
  expect_s3_class(res$sample_clusters, "data.frame")
  expect_equal(
    sort(unique(res$sample_clusters$cluster)),
    seq_len(4)
  )
  expect_equal(
    nrow(res$sample_clusters),
    length(unique(dat$Tumor_Sample_Barcode))
  )
  expect_equal(nrow(res$similarity_matrix), nrow(res$sample_clusters))
})

test_that("jaccard_cluster accepts a data frame input", {
  dat <- panel_subset()
  res <- jaccard_cluster(dat[, c("Tumor_Sample_Barcode", "Hugo_Symbol")],
                         n_clusters = 3, plot = FALSE)
  expect_equal(nrow(res$sample_clusters),
               length(unique(dat$Tumor_Sample_Barcode)))
})

test_that("MAF object input with unused factor levels clusters correctly", {
  # The merged MAF object stores its data as a data.table whose factor
  # columns keep unused levels (the full cohort); these must not leak into
  # table()/split() downstream (see .maf_to_dataframe()).
  sub <- panel_subset(20)
  expect_true(is.factor(LCBM_panel_maf@data$Tumor_Sample_Barcode))

  res <- jaccard_cluster(sub, n_clusters = 2, plot = FALSE)
  expect_equal(nrow(res$sample_clusters), 20)
  expect_equal(sort(unique(res$sample_clusters$cluster)), 1:2)

  genes <- cluster_top_genes(sub, res, top = 3)
  expect_true(all(c("cluster", "gene", "count") %in% colnames(genes)))
  expect_equal(nrow(genes), 6)
})

test_that("preprocess_maf_top_genes keeps the sample universe by default", {
  dat <- wes_subset()
  top <- preprocess_maf_top_genes(dat, max_features = 20)

  samples <- attr(top, "samples")
  expect_type(samples, "character")
  expect_setequal(samples, unique(as.character(dat$Tumor_Sample_Barcode)))
})

test_that("drop_empty_samples = TRUE restricts the sample universe", {
  dat <- wes_subset()
  top <- preprocess_maf_top_genes(dat, max_features = 20)
  strict <- preprocess_maf_top_genes(dat, max_features = 20,
                                     drop_empty_samples = TRUE)

  expect_equal(nrow(top), nrow(strict))
  expect_setequal(attr(strict, "samples"),
                  unique(as.character(strict$Tumor_Sample_Barcode)))
  expect_gte(length(attr(top, "samples")), length(attr(strict, "samples")))
})

# Small synthetic cohort used by the "empty sample" tests.
# Gene frequencies: A = 3, B = 2, C = 2, D = 1, Z = 1, so keeping the top 4
# genes drops Z and leaves S6 without any mutation in the retained genes
# (an all-zero similarity profile).
empty_sample_maf <- data.frame(
  Tumor_Sample_Barcode = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S5", "S6"),
  Hugo_Symbol          = c("A", "B", "A", "B", "A", "C", "C", "D", "Z"),
  stringsAsFactors     = FALSE
)

test_that("samples without mutation in the top genes keep an all-zero profile", {
  top <- preprocess_maf_top_genes(empty_sample_maf, max_features = 4)

  expect_setequal(attr(top, "samples"), c("S1", "S2", "S3", "S4", "S5", "S6"))
  expect_setequal(unique(top$Hugo_Symbol), c("A", "B", "C", "D"))
  expect_equal(length(unique(top$Tumor_Sample_Barcode)), 5L)

  sim <- jaccard_similarity_matrix(top)
  expect_equal(dim(sim), c(6L, 6L))
  expect_false(any(is.nan(sim)))

  # S6 has no mutation among the retained genes -> all-zero row and column.
  expect_true(all(sim["S6", ] == 0))
  expect_true(all(sim[, "S6"] == 0))

  # Retained legacy behaviour: two samples with identical gene sets score the
  # largest observed off-diagonal similarity (0.5 here), not 1.
  expect_equal(unname(sim["S1", "S2"]), 0.5^0.6)
  expect_equal(unname(sim["S1", "S3"]), (1/3)^0.6)
  expect_equal(unname(sim["S1", "S4"]), 0)

  strict <- preprocess_maf_top_genes(empty_sample_maf, max_features = 4,
                                     drop_empty_samples = TRUE)
  expect_equal(dim(jaccard_similarity_matrix(strict)), c(5L, 5L))
  expect_false("S6" %in% attr(strict, "samples"))
})

test_that("jaccard_cluster keeps empty samples by default and can drop them", {
  top <- preprocess_maf_top_genes(empty_sample_maf, max_features = 4)

  res <- jaccard_cluster(top, n_clusters = 2, clustering_method = "ward.D2",
                         plot = FALSE)
  expect_equal(nrow(res$sample_clusters), 6L)
  expect_equal(nrow(res$similarity_matrix), 6L)
  expect_true("S6" %in% res$sample_clusters$Tumor_Sample_Barcode)

  strict <- preprocess_maf_top_genes(empty_sample_maf, max_features = 4,
                                     drop_empty_samples = TRUE)
  res_strict <- jaccard_cluster(strict, n_clusters = 2,
                                clustering_method = "ward.D2", plot = FALSE)
  expect_equal(nrow(res_strict$sample_clusters), 5L)
  expect_false("S6" %in% res_strict$sample_clusters$Tumor_Sample_Barcode)
})

test_that("the shipped cohorts match the published analysis inputs", {
  expect_equal(length(unique(LCBM_panel_maf@data$Tumor_Sample_Barcode)), 358L)
  expect_equal(nrow(LCBM_panel_maf@data), 4691L)
  expect_equal(sort(unique(as.character(LCBM_panel_maf@data$cohort_name))),
               c("Huashan", "cbioportal"))

  expect_equal(length(unique(LCBM_WES_maf@data$Tumor_Sample_Barcode)), 157L)
  expect_equal(nrow(LCBM_WES_maf@data), 1585L)
  expect_equal(sort(unique(as.character(LCBM_WES_maf@data$cohort_name))),
               c("WES_SHANGHAI", "cellreport_WES"))

  # The annotation table covers every sample of both cohorts.
  samples <- c(as.character(LCBM_panel_maf@data$Tumor_Sample_Barcode),
               as.character(LCBM_WES_maf@data$Tumor_Sample_Barcode))
  expect_equal(sum(unique(samples) %in% LCBM_sample_annotation$Tumor_Sample_Barcode),
               length(unique(samples)))
})

test_that("the top-20 gene table matches the shipped expected results", {
  expected <- expected_csv("LCBM_top20_genes_panel.csv")
  counts <- sort(table(LCBM_panel_maf@data$Hugo_Symbol), decreasing = TRUE)

  expect_equal(expected$gene, names(counts)[seq_len(20)])
  expect_equal(expected$n_samples, as.integer(counts[seq_len(20)]))
  expect_equal(expected$gene[1], "TP53")
  expect_equal(expected$n_samples[1], 336L)
})

test_that("the shipped cluster tables hold the published 4-subtype split", {
  panel <- expected_csv("LCBM_cluster_panel.csv")
  wes   <- expected_csv("LCBM_cluster_WES.csv")

  expect_equal(nrow(panel), 342L)
  expect_equal(as.integer(table(panel$cluster)), c(178L, 65L, 73L, 26L))
  expect_equal(nrow(wes), 145L)
  expect_equal(as.integer(table(wes$cluster)), c(28L, 71L, 28L, 18L))

  # Empty samples are excluded, so a cluster table only covers a subset of the
  # MAF samples.
  expect_lt(nrow(panel), length(unique(LCBM_panel_maf@data$Tumor_Sample_Barcode)))
})
