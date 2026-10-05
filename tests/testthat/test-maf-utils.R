test_that("maf_common_columns returns shared columns", {
  common <- maf_common_columns(list(LCBM_panel_maf@data, LCBM_WES_maf@data))
  expect_true(all(c("Tumor_Sample_Barcode", "Hugo_Symbol") %in% common))
  # The two cohorts were harmonised, so they share far more than the two
  # mandatory columns.
  expect_true(all(c("VAF", "cohort_name", "DP") %in% common))
})

test_that("shape_maf harmonises and tags cohorts", {
  common <- maf_common_columns(list(LCBM_WES_maf@data))
  out <- shape_maf(LCBM_WES_maf@data, common, gene_use = c("TP53", "KRAS"),
                   cohort_name = "WES")
  expect_true(all(c("cohort_name") %in% colnames(out)))
  expect_true(all(out$Hugo_Symbol %in% c("TP53", "KRAS")))
  expect_true(all(out$cohort_name == "WES"))
})

test_that("plot_tmb_by_cluster returns a plot and wilcox table", {
  skip_if_not_installed("maftools")
  maf <- maftools::read.maf(wes_subset())
  res <- jaccard_cluster(wes_subset(), n_clusters = 3, plot = FALSE)
  out <- plot_tmb_by_cluster(maf, res)
  expect_s3_class(out$plot, "ggplot")
  expect_true(is.data.frame(out$wilcox))
})

test_that("plot_cluster_oncoplot subsets samples and draws", {
  skip_if_not_installed("maftools")
  maf <- maftools::read.maf(wes_subset())
  res <- jaccard_cluster(wes_subset(), n_clusters = 3, plot = FALSE)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  maf_sub <- plot_cluster_oncoplot(maf, res, cluster_id = 1, top = 5)
  expect_true(methods::is(maf_sub, "MAF"))
})
