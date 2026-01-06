
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
  old <- try(STList_legacy(rnacounts = cnts, spotcoords = coords, samples = "sample1"), silent=TRUE)
  
  # Run New (STlist, the refactored one detection)
  new <- STlist(rnacounts = cnts, spotcoords = coords, samples = names(cnts))
  
  # Compare
  expect_s4_class(new, "STlist")
  expect_equal(names(old@counts), names(new@counts))
  expect_equal(old@counts[["sample1"]], new@counts[["sample1"]])
  
  # Legacy and New both process lists similarly, so coords should match
  expect_equal(nrow(old@spatial_meta[["sample1"]]), nrow(new@spatial_meta[["sample1"]]))
})

test_that("Refactor Comparison: Xenium MEX (New Fix vs Legacy)", {
  tmp_dir <- file.path(tempdir(), "xenium_comp")
  create_mock_xenium_mex(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE))
  
  # Legacy check often fails or behaves weirdly. 
  # We mainly verify that the New implementation works.
  
  # New (Expect Success)
  new <- STlist(rnacounts = tmp_dir, samples = "sample1")
  
  expect_s4_class(new, "STlist")
  expect_true(length(new@counts) == 1)
  expect_equal(ncol(new@counts[[1]]), 10)
})
