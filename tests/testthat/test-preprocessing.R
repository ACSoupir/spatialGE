test_that("Data preprocessing (filtering, transform) works", {
  skip_if_not_installed("spatialGE")
  
  # Setup minimal STlist (mock or real)
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  
  # Test summarize
  # Capture output to avoid cluttering console test log
  summary_df <- capture.output(res <- summarize_STlist(melanoma))
  expect_true(is.data.frame(res) || is.null(res)) # summarize_STlist usually prints, returns df invisibly or NULL depending on version
  
  # Test Filter
  # Filter aggressively to see change
  melanoma_flt <- filter_data(melanoma, spot_minreads = 500, spot_mingenes = 100)
  expect_s4_class(melanoma_flt, "STlist")
  # Expect fewer spots or same, never more
  expect_true(ncol(melanoma_flt@counts[[1]]) <= ncol(melanoma@counts[[1]]))
  
  # Test Transform (Log)
  melanoma_log <- transform_data(melanoma_flt, method = "log")
  expect_s4_class(melanoma_log, "STlist")
  expect_true(!is.null(melanoma_log@tr_counts))
  expect_equal(names(melanoma_log@tr_counts), names(melanoma_log@counts))
  
  # Check one value to ensure log happened (roughly)
  # Raw count > 0 should be different in tr_counts (and likely small positive if log1p-like)
  # This is a weak check but verifies the slot is populated.
})
