test_that("Visualization functions return ggplot objects", {
  skip_if_not_installed("spatialGE")
  
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma, method='log')
  
  # STplot
  # Mocking gene presence: check checking top expressed genes first maybe?
  # Or just pick a gene we know exists. Melanoma usually expresses MLANA or typical Housekeepers.
  # Let's check first gene.
  gene_to_plot <- rownames(melanoma@tr_counts[[1]])[1]
  
  # Test gene plot
  p <- STplot(melanoma, genes=gene_to_plot)
  expect_true(is.list(p))
  expect_s3_class(p[[1]], "ggplot")
  
  # Test with user palette
  # "viridis" might not be in the khroma/RColorBrewer set supported by `color_parse` directly or needs a wrapper.
  # Let's use a known RColorBrewer palette like "Blues" or "Spectral"
  p2 <- STplot(melanoma, genes=gene_to_plot, color_pal="Blues")
  expect_s3_class(p2[[1]], "ggplot")
  
  # Pseudobulk Heatmap (requires pseudobulk first)
  # pseudobulk_samples requires >3 samples, our test set has 1. Skip for now.
  # melanoma <- pseudobulk_samples(melanoma)
  # ph <- pseudobulk_heatmap(melanoma)
  # expect_true(!is.null(ph))
})
