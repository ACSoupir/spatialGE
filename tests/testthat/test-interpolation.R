
test_that("Interpolation functions work", {
  # Load data
  data_dir <- file.path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1:2]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1:2]
  
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  
  # Test gene_interpolation
  # Interpolate top 2 genes for 1st sample
  st_obj <- gene_interpolation(st_obj, genes = 'top', top_n = 2, samples = 1, verbose = FALSE)
  
  expect_true(!is.null(st_obj@gene_krige))
  
  # Test STplot_interpolation
  # Need to know which genes were interpolated
  # They are stored in x@gene_krige keys
  interp_genes <- names(st_obj@gene_krige)
  
  p_interp <- STplot_interpolation(st_obj, genes = interp_genes[1], samples = 1)
  expect_type(p_interp, "list")
  expect_s3_class(p_interp[[1]], "ggplot")
})
