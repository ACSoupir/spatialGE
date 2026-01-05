
test_that("Enrichment and Gradient functions work", {
  # Load data
  data_dir <- file.path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  #clinical <- file.path(data_dir, "thrane_clinical.csv")

  clinical <- file.path(data_dir, "thrane_clinical.csv")
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  
  # --- STenrich Test ---
  # Mock gene sets
  # Get some gene names from first sample
  genes <- rownames(st_obj@tr_counts[[1]])
  gene_sets <- list(
    GS1 = genes[1:5],
    GS2 = genes[6:10]
  )
  
  # Run STenrich
  # Using domain=NULL to run on all spots
  # annot=NULL
  res_enrich <- STenrich(st_obj, gene_sets = gene_sets, samples = 1, reps = 10, verbose = FALSE)
  
  expect_type(res_enrich, "list")
  expect_true(length(res_enrich) > 0)
  expect_s3_class(res_enrich[[1]], "data.frame")
  
  # --- STgradient Test ---
  # Needs clustering
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = F)
  annot_col <- "stclust_spw0.05_k2"
  
  # Identify a reference cluster
  ref_cluster <- "1" # Assuming clusters are 1, 2...
  
  # Run STgradient
  # topgenes=50 for speed
  res_grad <- STgradient(st_obj, samples = 1, annot = annot_col, ref = ref_cluster, topgenes = 50, verbose = FALSE)
  
  expect_type(res_grad, "list")
  # It might be empty if no gradients found or no spots, but we expect a list
  
})
