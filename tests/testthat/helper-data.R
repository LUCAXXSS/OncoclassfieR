# pkgload::load_all() (used by devtools::test()) does not attach data/ to the
# search path, so the shipped datasets are always accessed explicitly.

# Real cohorts shipped with the package: these are what the tests exercise.
LCBM_panel_maf         <- OncoclassfieR::LCBM_panel_maf
LCBM_WES_maf           <- OncoclassfieR::LCBM_WES_maf
LCBM_sample_annotation <- OncoclassfieR::LCBM_sample_annotation

# The internal simulated datasets (`example_maf`, `example_survival`) are used
# by the function examples rather than by the tests; the edge cases that need a
# toy cohort define their own inline fixture (see `empty_sample_maf` in
# test-jaccard.R).

# Small slices of the real cohorts. Clustering 358 x 358 samples takes about
# half a minute, so the tests run the same code path on 60 panel / 30 WES
# samples, which keeps the suite well under a second per case while still
# going through the real data (factor barcodes, extra annotation columns).
cohort_subset <- function(maf, n) {
  dt <- as.data.frame(maf@data)
  keep <- unique(as.character(dt$Tumor_Sample_Barcode))[seq_len(n)]
  dt[as.character(dt$Tumor_Sample_Barcode) %in% keep, , drop = FALSE]
}

panel_subset <- function(n = 60) cohort_subset(LCBM_panel_maf, n)
wes_subset   <- function(n = 30) cohort_subset(LCBM_WES_maf, n)

# Survival table for the survival tests: the real clinical annotation,
# restricted to samples with an observed overall-survival time.
survival_data <- local({
  ann <- LCBM_sample_annotation
  ann <- ann[!is.na(ann$OS_time) & !is.na(ann$OS_status), , drop = FALSE]
  ann$OS_time   <- as.numeric(ann$OS_time)
  ann$OS_status <- as.integer(ann$OS_status)
  ann$cohort_name <- factor(ann$cohort_name)
  ann
})

# Expected analysis results shipped in inst/extdata/LCBM_expected.
expected_dir <- function() {
  system.file("extdata", "LCBM_expected", package = "OncoclassfieR")
}

expected_csv <- function(file) {
  utils::read.csv(file.path(expected_dir(), file), stringsAsFactors = FALSE)
}
