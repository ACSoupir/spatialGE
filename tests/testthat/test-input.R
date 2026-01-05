test_that("STlist creation from Thrane (Melanoma) data works", {
  skip_if_not_installed("spatialGE")
  
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")

  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  # Ensure we found files
  expect_true(length(count_files) > 0)
  expect_true(length(coord_files) > 0)
  expect_true(length(clin_file) > 0)
  
  # Create STlist using only first 2 samples to be faster
  melanoma <- STlist(
    rnacounts = count_files[1:2], 
    spotcoords = coord_files[1:2], 
    samples = clin_file
  )
  
  expect_s4_class(melanoma, "STlist")
  expect_equal(length(melanoma@counts), 2)
  expect_true(!is.null(melanoma@spatial_meta))
})

test_that("STlist creation handles missing inputs gracefully", {
    # Expect error or warning when inputs are wrong
    expect_error(STlist(rnacounts = "nonexistent.csv", spotcoords = "nonexistent.csv"))
})

# Note: Visium and SMI tests would go here, we can use the downloaded data if we implemented import_* tests.
# For now, focusing on the core STlist constructor which underpins imports.
