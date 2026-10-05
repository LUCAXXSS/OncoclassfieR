test_that("silhouette_analysis returns silhouette table and a plot", {
  res <- jaccard_cluster(panel_subset(), n_clusters = 4, plot = FALSE)
  out <- silhouette_analysis(res, k_range = 2:4)

  expect_type(out, "list")
  expect_named(out, c("silhouette", "plot"))
  expect_s3_class(out$silhouette, "data.frame")
  expect_named(out$silhouette, c("k", "silhouette"))
  expect_equal(nrow(out$silhouette), 3)
  expect_equal(out$silhouette$k, 2:4)
  expect_true(all(out$silhouette$silhouette >= -1 &
                  out$silhouette$silhouette <= 1))
  expect_s3_class(out$plot, "ggplot")
})

test_that("silhouette_analysis accepts a similarity matrix directly", {
  sim <- jaccard_similarity_matrix(wes_subset())
  out <- silhouette_analysis(sim, k_range = 2:3)

  expect_s3_class(out$silhouette, "data.frame")
  expect_equal(nrow(out$silhouette), 2)
  expect_equal(out$silhouette$k, 2:3)
  expect_s3_class(out$plot, "ggplot")
})
