test_that("Analysis functions (STclust, SThet) run without error", {
  skip_if_not_installed("spatialGE")
  
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # STclust
  # Use very small k or simple params to speed up
  # Note: STclust might be slow.
  melanoma <- STclust(melanoma, ks=2:3, ws=0.025)
  # Based on source, column should be stclust_spw0.025_k2 or k3
  expect_true(any(grepl("stclust_spw0.025_k", colnames(melanoma@spatial_meta[[1]]))))
  
  # SThet
  # Pick top varied gene to minimize time, or just 1 gene
  genes <- rownames(melanoma@tr_counts[[1]])[1:2]
  melanoma <- SThet(melanoma, genes=genes, method='moran')
  
  res <- get_gene_meta(melanoma, sthet_only=TRUE)
  expect_true(is.data.frame(res))
  expect_true(all(genes %in% res$gene))
  
  # compare_SThet
  # Should return a plot (ggplot/ggarrange)
  p_cmp <- compare_SThet(melanoma, genes=genes)
  expect_true(is(p_cmp, "ggplot") || is(p_cmp, "list"))
  
})
