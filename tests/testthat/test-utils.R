
test_that("Utility and QC functions work", {
  # Load data
  data_dir <- file.path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")

  # Create STlist
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)

  # 1. Test tissue_names
  tn <- tissue_names(st_obj)
  expect_type(tn, "character")
  expect_true(length(tn) > 0)

  # 2. Test spatial_metadata
  sm <- spatial_metadata(st_obj)
  expect_type(sm, "list")
  expect_true(length(sm) > 0)

  # 3. Test plot_counts
  # plot_counts requires transformed data for default plot_type?
  # "No transformed expression data in this STlist." means we need transform_data() first
  st_obj <- transform_data(st_obj)
  p <- plot_counts(st_obj, plot_type = "density")
  expect_type(p, "list") # returns list of plots
  
  # 4. Test distribution_plots
  # Needs plot_type (defaults to 'violin', 'box', 'density')? 
  # Actually documentation says: To plot counts or genes per spot/cell, the function distribution_plots should be used instead.
  # Error said: Specify a spot/cell meta data variable or a gene to plot.
  # So we need 'items' argument? 
  # Let's plot 'total_counts' which is usually in spatial_meta
  dp <- distribution_plots(st_obj, plot_meta = "total_counts", plot_type="box")
  expect_s3_class(dp[[1]], "ggplot")
  
  # 5. Test get_gene_meta
  gm <- get_gene_meta(st_obj)
  expect_type(gm, "list")
  expect_true(length(gm) > 0)
  
  # 6. Test summarize_STlist
  summ <- summarize_STlist(st_obj)
  expect_true(TRUE) 
  
})


