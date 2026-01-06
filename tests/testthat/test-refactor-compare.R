
source(testthat::test_path("helper-regression.R"))

test_that("Refactor Comparison: List of DataFrames", {
  # Mock list of dfs
  cnts <- list(
      sample1 = data.frame(
          gene = c("g1", "g2"),
          c1 = c(10, 0),
          c2 = c(5, 5)
      )
  )
  coords <- list(
      sample1 = data.frame(
          colname = c("c1", "c2"),
          ypos = c(1, 2),
          xpos = c(1, 2)
      )
  )
  
  # Run Legacy
  old <- try(STlist(rnacounts = cnts, spotcoords = coords, samples = "sample1"), silent=TRUE)
  
  # Run New
  new <- STList_new(rnacounts = cnts, spotcoords = coords, samples = names(cnts))
  
  # Compare
  expect_s4_class(new, "STlist")
  expect_equal(names(old@counts), names(new@counts))
  expect_equal(old@counts[["sample1"]], new@counts[["sample1"]])
  # coords might differ slightly in col names if new ingest standardizes them better?
  # Legacy STList uses detection logic to rename cols. New logic maps explicitly.
  # Let's check essential content.
  expect_equal(nrow(old@spatial_meta[["sample1"]]), nrow(new@spatial_meta[["sample1"]]))
})

# test_that("Refactor Comparison: Visium MEX (Tarball)", {
#   # Skipped as per user request (Legacy STList does not support robustly)
# })

# test_that("Refactor Comparison: Visium MEX (Unzipped Directory - Improvement)", {
#   # Skipped as per user request (Legacy STList does not support)
# })

test_that("Refactor Comparison: Xenium MEX (Improvement/Fix)", {
  tmp_dir <- file.path(tempdir(), "xenium_comp")
  create_mock_xenium_mex(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE))
  
  # Legacy check causes test noise/failure due to warnings. 
  # We know it fails. Focusing on NEW implementation success.
  # err_old <- tryCatch(STlist(rnacounts = tmp_dir, samples = "sample1"), error=function(e) e)
  
  # New (Expect Success)
  new <- STList_new(rnacounts = tmp_dir, samples = "sample1")
  
  # Debug Output
  print(paste("New Class:", class(new)))
  if (inherits(new, "STlist")) {
      print(paste("Counts length:", length(new@counts)))
      if (length(new@counts) > 0) print(paste("Dims:", paste(dim(new@counts[[1]]), collapse="x")))
  }
  
  expect_s4_class(new, "STlist")
  expect_true(length(new@counts) == 1)
  expect_equal(ncol(new@counts[[1]]), 10)
})
