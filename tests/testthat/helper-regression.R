
# Helper to create mock Visium MEX directory
create_mock_visium_mex <- function(dir_path, sample_name="sample1") {
  # Dir structure: sample/spatial/tissue_positions_list.csv (or .parquet)
  #                sample/filtered_feature_bc_matrix/matrix.mtx.gz
  #                sample/filtered_feature_bc_matrix/barcodes.tsv.gz
  #                sample/filtered_feature_bc_matrix/features.tsv.gz
  
  if(dir.exists(dir_path)) unlink(dir_path, recursive=TRUE)
  dir.create(dir_path, recursive=TRUE)
  
  # 1. Create spatial coords
  spatial_dir <- file.path(dir_path, "spatial")
  dir.create(spatial_dir)
  
  # tissue_positions_list.csv: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
  barcodes <- paste0("BC", 1:10)
  coords <- data.frame(
    barcode = barcodes,
    in_tissue = 1,
    array_row = 1:10,
    array_col = 1:10,
    pxl_row = (1:10)*10,
    pxl_col = (1:10)*10
  )
  write.table(coords, file.path(spatial_dir, "tissue_positions_list.csv"), 
              sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # 2. Create Matrices
  mat_dir <- file.path(dir_path, "filtered_feature_bc_matrix")
  dir.create(mat_dir)
  
  # Barcodes
  gz_bc <- gzfile(file.path(mat_dir, "barcodes.tsv.gz"), "w")
  writeLines(barcodes, gz_bc)
  close(gz_bc)
  
  # Features: id, name, type
  genes <- paste0("Gene", 1:5)
  feats <- data.frame(
      id = genes,
      name = genes,
      type = "Gene Expression"
  )
  gz_ft <- gzfile(file.path(mat_dir, "features.tsv.gz"), "w")
  write.table(feats, gz_ft, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(gz_ft)
  
  # Matrix (.mtx)
  # Header
  header <- c("%%MatrixMarket matrix coordinate integer general", "%", 
              paste(length(genes), length(barcodes), length(genes)*length(barcodes)))
  
  # Data: row col val
  # All 1s
  rows <- rep(1:length(genes), times=length(barcodes))
  cols <- rep(1:length(barcodes), each=length(genes))
  vals <- rep(1, length(rows))
  vals <- rep(1, length(rows))
  df_mat <- data.frame(rows, cols, vals)
  
  gz_mx <- gzfile(file.path(mat_dir, "matrix.mtx.gz"), "w")
  writeLines(header, gz_mx)
  write.table(df_mat, gz_mx, sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(gz_mx)
  
  # Tar the matrix directory (STList expects .tar.gz for MEX)
  # tar command available on mac/linux
  current_wd <- getwd()
  setwd(dir_path)
  
  # Visium: filtered_feature_bc_matrix.tar.gz (Strict requirement for detection)
  utils::tar("filtered_feature_bc_matrix.tar.gz", "filtered_feature_bc_matrix", compression="gzip")
  unlink("filtered_feature_bc_matrix", recursive=TRUE)
  
  setwd(current_wd)
  
  return(dir_path)
}

# Helper to create mock Xenium MEX directory
create_mock_xenium_mex <- function(dir_path, sample_name="sample1") {
  # Dir structure: sample/cell_feature_matrix.tar.gz
  #                sample/cells.csv (or parquet)
  
  if(dir.exists(dir_path)) unlink(dir_path, recursive=TRUE)
  dir.create(dir_path, recursive=TRUE)
  
  # 1. Coordinate file
  barcodes <- paste0("Cell", 1:10)
  coords <- data.frame(
      cell_id = barcodes,
      x_centroid = (1:10)*10,
      y_centroid = (1:10)*10,
      transcript_counts = 100,
      control_probe_counts = 5,
      control_codeword_counts = 5,
      total_counts = 110
  )
  write.csv(coords, file.path(dir_path, "cells.csv"), row.names=FALSE)
  
  # 2. Matrix
  mat_dir <- file.path(dir_path, "cell_feature_matrix")
  dir.create(mat_dir)
  
  gz_bc <- gzfile(file.path(mat_dir, "barcodes.tsv.gz"), "w")
  writeLines(barcodes, gz_bc)
  close(gz_bc)
  
  genes <- paste0("Gene", 1:5)
  feats <- data.frame(id=genes, name=genes, type="Gene Expression")
  gz_ft <- gzfile(file.path(mat_dir, "features.tsv.gz"), "w")
  write.table(feats, gz_ft, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(gz_ft)
  
  header <- c("%%MatrixMarket matrix coordinate integer general", "%", 
              paste(length(genes), length(barcodes), length(genes)*length(barcodes)))
  rows <- rep(1:length(genes), times=length(barcodes))
  cols <- rep(1:length(barcodes), each=length(genes))
  vals <- rep(1, length(rows))
  df_mat <- data.frame(rows, cols, vals)
  gz_mx <- gzfile(file.path(mat_dir, "matrix.mtx.gz"), "w")
  writeLines(header, gz_mx)
  write.table(df_mat, gz_mx, sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(gz_mx)
  
  # Tar it
  current_wd <- getwd()
  setwd(dir_path)
  # cell_feature_matrix.tar.gz
  utils::tar("cell_feature_matrix.tar.gz", "cell_feature_matrix", compression="gzip")
  unlink("cell_feature_matrix", recursive=TRUE)
  setwd(current_wd)

  return(dir_path)
}
