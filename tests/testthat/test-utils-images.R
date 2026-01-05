
test_that("load_images and plot_image work with real data", {
  # Use TNBC Bassiouni dataset if available (it has images)
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")

  # Load one sample from TNBC
  # sample_092a
  sample_dir <- file.path(data_dir, "sample_092a")
  if(dir.exists(sample_dir)){
      # Create STlist
      # For Visium, rnacounts is the dir
      st_obj <- STlist(rnacounts = sample_dir, samples = "sample_092a")
      
      # Image path
      img_path <- file.path(sample_dir, "spatial", "GSM6433585_092A_tissue_hires_image.png")
      skip_if_not(file.exists(img_path), "Image file not found")
      
      # 1. Test load_images
      st_obj <- load_images(st_obj, images = img_path)
      expect_true(!is.null(st_obj@misc[['sp_images']][['sample_092a']]))
      
      # 2. Test plot_image
      p_img <- plot_image(st_obj)
      expect_type(p_img, "list")
      expect_s3_class(p_img[[1]], "ggplot")
  }
})

test_that("load_images handles missing images gracefully", {
    # Keep the negative test with melanoma data
    data_dir <- file.path("data", "melanoma_thrane")
    counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
    coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
    clinical <- file.path(data_dir, "thrane_clinical.csv")
    st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
    
    # It seems to not warn or error, just return object or print message
    res_img <- load_images(st_obj, images = "dummy.jpg")
    expect_true(is(res_img, "STlist"))
    expect_error(plot_image(st_obj), "No tissue images available")
})
