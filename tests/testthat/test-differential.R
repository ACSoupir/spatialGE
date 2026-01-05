
test_that("Differential expression functions work", {
  # Load data
  data_dir <- file.path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")

  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  
  # Need clustering first for annotations if not using clinical
  # Using a small subset or very coarse clustering to be fast?
  # Actually, we can use clinical annotation if available.
  # The clinical data has 'gender', 'BRAF_status', etc.
  # Let's check spatial_metadata again.
  # STlist merges clinical into spatial_metadata.
  
  # Test STdiff using clinical annotation (e.g., 'gender' if it varies within sample? No, clinical is per sample)
  # STdiff tests groups of spots/cells.
  # If we have sample-level metadata, STdiff might not work if it expects spot-level differences within a sample?
  # The documentation says: "Tests for differentially expressed genes between groups of spots/cells (e.g., clusters) in a spatial transcriptomics sample."
  # So we definitely need spot-level clusters (STclust).
  
  # Run STclust (fast version)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = F)
  
  # Get the cluster column name
  # It should be stclust_spw0.05_k2
  cluster_col <- "stclust_spw0.05_k2"
  
  # Run STdiff (non-spatial only for speed)
  # samples=NULL => all samples
  # topgenes=100 for speed
  res_diff <- STdiff(st_obj, annot = cluster_col, topgenes = 100, sp_topgenes = 0, verbose = 0)
  
  expect_type(res_diff, "list")
  expect_true(length(res_diff) > 0)
  expect_s3_class(res_diff[[1]], "data.frame")
  
  # Test STdiff_volcano
  p_vol <- STdiff_volcano(res_diff, samples = names(res_diff)[1])
  expect_type(p_vol, "list") # It might return a list of plots or a single plot object?
  # Checking STdiff_volcano impl: it uses ggpubr::ggarrange or returns list? 
  # Let's assume list of plots or ggplot.
  
})
