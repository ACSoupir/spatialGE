# Helper script to download test data if not present
# This script is automatically sourced by testthat before running tests

setup_test_data <- function() {
  data_dir <- testthat::test_path("data")
  
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Define datasets to download
  datasets <- list(
    tnbc = list(
      url = "https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/tnbc_bassiouni.zip",
      folder_name = "tnbc_bassiouni"
    ),
    melanoma = list(
      url = "https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip",
      folder_name = "melanoma_thrane"
    ),
    nsclc = list(
      url = "https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/nanostring_nsclc.zip",
      folder_name = "nanostring_nsclc"
    )
  )
  
  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    target_dir <- file.path(data_dir, ds$folder_name)
    
    # Check if data already exists
    if (!dir.exists(target_dir)) {
      message(sprintf("Downloading %s data...", ds_name))
      
      zip_file <- file.path(data_dir, paste0(ds$folder_name, ".zip"))
      
      tryCatch({
        download.file(ds$url, destfile = zip_file, mode = "wb", quiet = TRUE)
        unzip(zip_file, exdir = data_dir)
        unlink(zip_file) # Clean up zip file
        message(sprintf("Successfully setup %s data.", ds_name))
      }, error = function(e) {
        warning(sprintf("Failed to download or unzip %s data: %s", ds_name, e$message))
      })
    }
  }
}

# Run setup
setup_test_data()
