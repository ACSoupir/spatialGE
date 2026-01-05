
test_that("Pseudobulk functions work with real data", {
  # Use TNBC Bassiouni dataset (Visium) which has > 4 samples
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")

  # Setup paths for Visium data (folders containing .h5)
  # The helper script unzips them into sample_XXX folders
  # But STlist for Visium expects paths to directories containing the .h5 files
  
  # List all sample directories
  sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  # Filter only those starting with "sample_"
  sample_dirs <- grep("sample_", sample_dirs, value = TRUE)
  
  # Ensure we have at least 4
  expect_true(length(sample_dirs) >= 4)
  
  # Clinical data
  clinical <- file.path(data_dir, "bassiouni_clinical.csv")
  
  # STlist for Visium: rnacounts = vector of paths to directories
  # output structure: samples are usually named optionally or inferred
  # We can provide a samples vector or file
  
  # Create STlist
  # Visium data loading might be slower, so we use a subset if needed, but we need 4 samples.
  use_samples <- sample_dirs[1:4]
  
  # To avoid warnings about missing samples (because clinical file has 8 samples but we load 4),
  # we create a subset clinical file.
  clin_df <- read.csv(clinical)
  
  # Extract sample IDs from folder names (basename)
  # The folder names are "sample_092a", etc.
  use_ids <- basename(use_samples)
  
  # Check if IDs match clinical
  # The clinical file has 'sample_id' column
  if("sample_id" %in% colnames(clin_df)) {
      clin_sub <- clin_df[clin_df$sample_id %in% use_ids, ]
      
      tmp_clin <- tempfile(fileext = ".csv")
      write.csv(clin_sub, tmp_clin, row.names = FALSE)
      
      # Use subset clinical
      st_obj <- STlist(rnacounts = use_samples, samples = tmp_clin)
  } else {
      # Fallback if column names differ
      st_obj <- STlist(rnacounts = use_samples, samples = clinical)
  }
  
  # Verify we loaded enough samples
  expect_true(length(st_obj@counts) >= 4)
  
  # Test pseudobulk_samples
  # Use small max_var_genes to avoid errors if few genes match or for speed
  st_obj <- pseudobulk_samples(st_obj, max_var_genes = 500)
  
  expect_true(!is.null(st_obj@misc[['scaled_pbulk_mtx']]))
  
  # Test pseudobulk_dim_plot
  p_pca <- pseudobulk_dim_plot(st_obj, dim='pca')
  expect_s3_class(p_pca, "ggplot")
  
  # Test pseudobulk_heatmap
  # This was the original goal of the user (to fix the skipped test)
  p_hm <- pseudobulk_heatmap(st_obj, hm_display_genes = 20)
  # Heatmap returns ComplexHeatmap object or similar (drawn)
  expect_true(is(p_hm, "Heatmap") || is(p_hm, "HeatmapList"))
})
