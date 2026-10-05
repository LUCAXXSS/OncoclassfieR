test_that("cluster_top_genes ranks genes per cluster", {
  dat <- panel_subset()
  res <- jaccard_cluster(dat, n_clusters = 4, plot = FALSE)
  genes <- cluster_top_genes(dat, res, top = 5)

  expect_s3_class(genes, "data.frame")
  expect_named(genes, c("cluster", "gene", "count"))
  expect_equal(nrow(genes), 4 * 5)
  expect_true(all(genes$count >= 1))
})

test_that("cluster_top_genes accepts a plain cluster data frame", {
  dat <- wes_subset()
  res <- jaccard_cluster(dat, n_clusters = 3, plot = FALSE)
  genes <- cluster_top_genes(dat,
                             res$sample_clusters, top = 3)
  expect_equal(nrow(genes), 3 * 3)
})

test_that("plot_cluster_piechart returns a ggplot", {
  dat <- wes_subset()
  res <- jaccard_cluster(dat, n_clusters = 3, plot = FALSE)
  genes <- cluster_top_genes(dat, res, top = 3)
  p <- plot_cluster_piechart(genes, top = 3)
  expect_s3_class(p, "ggplot")
})

test_that("pie charts retain every slice at the 100-percent boundary", {
  genes <- expected_csv("LCBM_top_genes_WES.csv")
  p <- plot_cluster_piechart(genes, top = 3)
  built <- ggplot2::ggplot_build(p)
  bars <- built$data[[1]]
  labels <- built$data[[2]]

  expect_equal(nrow(bars), nrow(p$data))
  expect_true(all(is.finite(bars$ymin)))
  expect_true(all(is.finite(bars$ymax)))
  expect_true(all(is.finite(labels$y)))
  expect_true(all(bars$ymax > bars$ymin))
  expect_equal(as.numeric(tapply(bars$ymax, bars$PANEL, max)), rep(100, 4))
  expect_no_warning(ggplot2::ggplotGrob(p))
})
