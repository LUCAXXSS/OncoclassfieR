test_that("fit_km returns a survfit object", {
  fit <- fit_km(survival_data, time = "OS_time", event = "OS_status",
                group = "cohort_name")
  expect_s3_class(fit, "survfit")
})

test_that("fit_cox returns a coxph object", {
  fit <- fit_cox(survival_data, time = "OS_time", event = "OS_status",
                 group = "cohort_name")
  expect_s3_class(fit, "coxph")
})

test_that("plot_km_curves returns a ggsurvplot", {
  fit <- fit_km(survival_data, time = "OS_time", event = "OS_status",
                group = "cohort_name")
  p <- plot_km_curves(fit, survival_data)
  expect_true(inherits(p, "ggsurvplot"))
})

test_that("missing columns are detected", {
  expect_error(fit_km(survival_data, time = "OS_time", event = "OS_status",
                      group = "missing_col"),
               "missing column")
})
