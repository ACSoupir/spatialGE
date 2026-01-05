
test_that("Pseudobulk functions work", {
  # Mock 4 samples to satisfy requirement
  data_dir <- file.path("data", "melanoma_thrane")
  
  # existing files
  c1 <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  m1 <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  
  # Create temp dir for mock data
  mock_dir <- tempfile()
  dir.create(mock_dir)
  on.exit(unlink(mock_dir, recursive = TRUE))
  
  # Create 4 dummy samples
  samp_names <- paste0("Sample", 1:4)
  c_files <- file.path(mock_dir, paste0(samp_names, "_counts.tsv"))
  m_files <- file.path(mock_dir, paste0(samp_names, "_mapping.tsv"))
  
  for(i in 1:4){
    # Read template counts
    # The counts file has genes as rows, spots as columns (usually)
    # The header might have "Gene" or empty string for the first column
    # uses check.names=FALSE to prevent X prefixing if spots start with numbers
    df <- read.table(c1, header=TRUE, row.names=NULL, sep="\t", fill=TRUE, check.names=FALSE)
    # Ensure first column is genes (unique)
    if(ncol(df) > 1) {
       df <- df[!duplicated(df[,1]),]
       rownames(df) <- df[,1]
       df <- df[,-1]
       
       # Save original column names
       orig_cols <- colnames(df)
       
       # Convert to numeric to add noise safely
       df[] <- lapply(df, function(x) as.numeric(as.character(x)))
       df[is.na(df)] <- 0
       
       # Add some noise to ensure variance across samples
       set.seed(i)
       noise <- round(rnorm(nrow(df) * ncol(df), mean=0, sd=1))
       df <- df + noise
       df[df < 0] <- 0 # Ensure non-negative
       
       # Write modified counts
       # STlist expects genes as rows. 
       # read.table with row.names=NULL read the gene names as first column.
       # we set them as rownames.
       
       # When writing, we need to match the input format expected by STlist
       # If original file had gene names in first column (header-less or named), we need to replicate
       # The error said "ROI/spot/cell IDs ... do not match".
       # This implies the columns of the count matrix must match the rows of the coordinate matrix.
       # The coordinate matrix (m_files) is just copied from m1.
       # So colnames(df) must match the spot IDs in m1.
       
       # We stripped the first column (genes), so the remaining columns should be spots.
       # We must ensure they match m1 spot IDs.
       
       # Let's peek at m1 to see what format it expects?
       # But simpler: we just need to preserve the column names from original read exactly.
       # 'orig_cols' contains the gene column name (index 1) and spot columns.
       # We removed index 1.
       
       # Write with tab sep
       write.table(df, c_files[i], sep="\t", quote=FALSE, col.names=NA)
    } else {
       # Fallback if file read failed oddly, just copy
       file.copy(c1, c_files[i])
    }
    
    file.copy(m1, m_files[i])
  }
  
  # Create STlist
  st_obj <- STlist(rnacounts = c_files, spotcoords = m_files, samples = samp_names)
  
  # Test pseudobulk_samples
  # Use small max_var_genes to avoid errors if few genes
  st_obj <- pseudobulk_samples(st_obj, max_var_genes = 100)
  
  expect_true(!is.null(st_obj@misc[['scaled_pbulk_mtx']]))
  
  # Test pseudobulk_dim_plot
  p_pca <- pseudobulk_dim_plot(st_obj, dim='pca')
  expect_s3_class(p_pca, "ggplot")
  
  # Test pseudobulk_heatmap
  p_hm <- pseudobulk_heatmap(st_obj, hm_display_genes = 10)
  # Heatmap returns ComplexHeatmap object or similar (drawn)
  expect_true(is(p_hm, "Heatmap") || is(p_hm, "HeatmapList"))
})
