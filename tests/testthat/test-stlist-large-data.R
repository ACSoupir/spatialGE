
test_that("STlist Ingests Large Seurat Data (Soupir et al)", {
  seurat_path <- testthat::test_path("data/soupir_etal_seurat/seurat_object.Rds")
  skip_if_not(file.exists(seurat_path), "Large Seurat data not found")
  
  # Load Seurat Object
  seurat_obj <- readRDS(seurat_path)
  
  # Ingest
  # STlist detection should handle Seurat object passed as rnacounts
  st <- STlist(rnacounts=seurat_obj, verbose=FALSE)
  
  expect_s4_class(st, "STlist")
  # Expect 3 samples as seen in verification (Run5448.RCC3, RCC4, RCC5)
  expect_equal(length(st@counts), 3)
  expect_true(all(c("Run5448.RCC3", "Run5448.RCC4", "Run5448.RCC5") %in% names(st@counts)))
})

test_that("STlist Ingests Large Xenium Data (Pancreas)", {
  xenium_path <- testthat::test_path("data/xenium_pancreas/Xenium_V1_human_Pancreas_FFPE_outs")
  skip_if_not(file.exists(xenium_path), "Large Xenium data not found")
  
  # Ingest
  # Pass directory path
  st <- STlist(rnacounts=xenium_path, samples="xenium_pancreas", verbose=FALSE)
  
  expect_s4_class(st, "STlist")
  expect_equal(length(st@counts), 1)
  expect_equal(names(st@counts), "xenium_pancreas")
  
  # Check dimensions (verification showed ~100k spots, 300 genes)
  # Just checking not empty
  expect_gt(nrow(st@counts[[1]]), 0)
  expect_gt(ncol(st@counts[[1]]), 0)
})
