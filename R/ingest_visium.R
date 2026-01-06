#' Ingest Visium Data
#' @param source InputSource object of type 'visium'
#' @export
dispatch_ingest.source_visium <- function(source) {
  path <- source$rna
  
  # Determine format if not explicit (directory structure vs specific file?)
  # For Visium source, path is usually a directory.
  
  # 1. Look for H5 (raw or filtered)
  h5_files <- list.files(path, pattern = "filtered_feature_bc_matrix.h5$", full.names = TRUE, recursive = TRUE)
  if (length(h5_files) == 0) {
      h5_files <- list.files(path, pattern = "raw_feature_bc_matrix.h5$", full.names = TRUE, recursive = TRUE)
  }
  
  if (length(h5_files) > 0) {
      # Prefer H5
      return(ingest_visium_h5(h5_files[1], source$coords, source$samples))
  }
  
  # 2. Look for MEX (filtered or raw)
  # Can be folder or tar.gz
  mex_dirs <- list.files(path, pattern = "filtered_feature_bc_matrix$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  if (length(mex_dirs) == 0) {
      mex_dirs <- list.files(path, pattern = "raw_feature_bc_matrix$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  }
  
  # Also check for tar.gz (legacy STList requirement)
  mex_tars <- list.files(path, pattern = "filtered_feature_bc_matrix.tar.gz$", full.names = TRUE, recursive = TRUE)
  
  if (length(mex_dirs) > 0) {
      return(ingest_visium_mex_dir(mex_dirs[1], source$coords, source$samples))
  } else if (length(mex_tars) > 0) {
      # Extract tar to temp loc
      tmp_ex <- file.path(tempdir(), paste0("vis_ex_", basename(path)))
      utils::untar(mex_tars[1], exdir=tmp_ex)
      # Find the dir inside
      inner_dirs <- list.files(tmp_ex, full.names=TRUE, include.dirs=TRUE)
      # Usually simple extraction creates the dir
      # Pass path as original_root
      return(ingest_visium_mex_dir(inner_dirs[1], source$coords, source$samples, original_root=path))
  }
  
  stop("Could not find valid Visium H5 or MEX data in: ", path)
}

#' Ingest Generic H5 Details (Visium or Xenium or 10x)
#' @param source InputSource object of type 'h5_10x'
#' @export
dispatch_ingest.source_h5_10x <- function(source) {
  path <- source$rna
  # H5 file path detected
  
  # Heuristic 1: Check companion files (coordinates) to distinguish
  # If coords provided, check extension or content
  is_xenium <- FALSE
  is_visium <- FALSE
  
  cpath <- source$coords
  if(!is.null(cpath)) {
      if(grepl("cells\\.parquet$|cells\\.csv(\\.gz)?$", cpath)) is_xenium <- TRUE
      if(grepl("tissue_positions.*\\.csv", cpath)) is_visium <- TRUE
  }
  
  if(!is_xenium && !is_visium) {
     # Check siblings
     parent <- dirname(path)
     if(length(list.files(parent, pattern="cells\\.parquet|cells\\.csv")) > 0) is_xenium <- TRUE
     if(length(list.files(parent, pattern="tissue_positions.*\\.csv")) > 0) is_visium <- TRUE
  }
  
  # If still ambiguous, check filename
  if(!is_xenium && !is_visium) {
      if(grepl("cell_feature_matrix\\.h5$", path)) is_xenium <- TRUE
      if(grepl("filtered_feature_bc_matrix\\.h5$|raw_feature_bc_matrix\\.h5$", path)) is_visium <- TRUE
  }
  
  # Dispatch
  # Note: ingest_xenium_h5 is available in package namespace (defined in ingest_xenium.R)
  if(is_xenium) {
      return(ingest_xenium_h5(path, source$coords, source$samples))
  } else {
      # Default to Visium (standard 10x H5)
      return(ingest_visium_h5(path, source$coords, source$samples))
  }
}

#' Internal H5 Reader
ingest_visium_h5 <- function(h5_path, coords_path=NULL, sample_name=NULL) {
  # Logic adapted from import_visium_h5
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("hdf5r required")
  
  h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(h5_file$close_all())
  
  # 10x standard: matrix/data, matrix/indices, etc.
  # Or root group
  root <- names(h5_file)[1]
  grp <- h5_file[[root]]
  
  # Read sparse matrix components
  data <- grp[["data"]][]
  indices <- grp[["indices"]][] 
  indptr <- grp[["indptr"]][]
  shape <- grp[["shape"]][]
  
  # Load features/barcodes
  if("features" %in% names(grp)) {
      feat_grp <- grp[["features"]]
      genes <- feat_grp[["name"]][]
      ids <- feat_grp[["id"]][]
  } else if ("genes" %in% names(grp)) { # Older format?
      genes <- grp[["genes"]][]
      ids <- genes # fallback
  } else {
      stop("Could not find features/genes in H5")
  }
  
  barcodes <- grp[["barcodes"]][]
  
  # Construct Matrix
  counts <- Matrix::sparseMatrix(
      x = as.numeric(data),
      i = indices,
      p = indptr,
      dims = shape,
      index1 = FALSE
  )
  rownames(counts) <- make.unique(genes)
  colnames(counts) <- barcodes
  
  # Load Coordinates
  # If coords_path is NULL, look relative to h5 (usually ../spatial/tissue_positions_list.csv)
  if (is.null(coords_path)) {
      parent <- dirname(dirname(h5_path)) # ../../ relative to h5 if nested?
      # Usually structure: sample/outs/filtered...h5
      #                    sample/outs/spatial/tissue...csv
      # Or: sample/filtered...h5
      #     sample/spatial/...
      
      # Try looking nearby
      candidates <- c(
          file.path(dirname(h5_path), "spatial", "tissue_positions_list.csv"),
          file.path(dirname(dirname(h5_path)), "spatial", "tissue_positions_list.csv")
      )
      
      # Also check for non-standard names (e.g. GSM prefix) in spatial dirs
      # Check strictly inside 'spatial' folders if they exist
      spatial_dir_1 <- file.path(dirname(h5_path), "spatial")
      if(dir.exists(spatial_dir_1)) {
        candidates <- c(candidates, list.files(spatial_dir_1, pattern="tissue_positions.*\\.csv", full.names=TRUE))
      }
      
      spatial_dir_2 <- file.path(dirname(dirname(h5_path)), "spatial")
      if(dir.exists(spatial_dir_2)) {
        candidates <- c(candidates, list.files(spatial_dir_2, pattern="tissue_positions.*\\.csv", full.names=TRUE))
      }
      
      valid <- candidates[file.exists(candidates)]
      if(length(valid) > 0) coords_path <- valid[1]
  }
  
  coords_df <- NULL
  if (!is.null(coords_path) && file.exists(coords_path)) {
      # Visium standard: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
      # header=FALSE usually
      raw_co <- data.table::fread(coords_path, header=FALSE)
      # Auto-detect if header exists (check if col1 is 'barcode')
      if (raw_co[1,1] == "barcode") {
          colnames(raw_co) <- as.character(raw_co[1,])
          raw_co <- raw_co[-1,]
      } else {
          colnames(raw_co) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
      }
      
      # Filter in_tissue
      raw_co <- raw_co[raw_co$in_tissue == 1, ]
      
      # Clean for STList (barcode, ypos, xpos) -> (barcode, imagerow, imagecol)
      coords_df <- data.frame(
          barcode = raw_co$barcode,
          imagerow = as.numeric(raw_co$pxl_row),
          imagecol = as.numeric(raw_co$pxl_col)
      )
      
      # Align
      common <- intersect(colnames(counts), coords_df$barcode)
      counts <- counts[, common, drop=FALSE]
      coords_df <- coords_df[match(common, coords_df$barcode), ]
  }
  
  return(list(counts=counts, coords=coords_df))
}

#' Internal MEX Directory Reader
ingest_visium_mex_dir <- function(mex_dir, coords_path=NULL, sample_name=NULL, original_root=NULL) {
  # Matrix::readMM for matrix.mtx
  mtx_path <- file.path(mex_dir, "matrix.mtx")
  if (!file.exists(mtx_path)) {
      mtx_path <- paste0(mtx_path, ".gz")
  }
  
  if (!file.exists(mtx_path)) {
      files_in_dir <- paste(list.files(mex_dir), collapse=", ")
      stop("Matrix file not found: ", mex_dir, "\nExpected matrix.mtx or matrix.mtx.gz.\nFiles found: ", files_in_dir)
  }
  
  counts <- read_mtx_safe(mtx_path)

  
  # Features/Genes
  feat_path <- file.path(mex_dir, "features.tsv")
  if (!file.exists(feat_path)) feat_path <- paste0(feat_path, ".gz")
  # Try genes.tsv if features not found
  if (!file.exists(feat_path)) {
       feat_path <- file.path(mex_dir, "genes.tsv")
       if (!file.exists(feat_path)) feat_path <- paste0(feat_path, ".gz")
  }
  
  feats <- data.table::fread(feat_path, header=FALSE)
  # Col 2 is typically gene name associated with rows
  genes <- feats[[2]]
  rownames(counts) <- make.unique(genes)
  
  # Barcodes
  bc_path <- file.path(mex_dir, "barcodes.tsv")
  if (!file.exists(bc_path)) bc_path <- paste0(bc_path, ".gz")
  bcs <- data.table::fread(bc_path, header=FALSE)[[1]]
  colnames(counts) <- bcs
  
  # Coordinates (Logic tailored to finding them relative to the input source root,
  # but here we are deep in MEX dir. Strategy: Look up 2 levels usually?)
  if (is.null(coords_path)) {
      # Heuristic search up tree
      root <- dirname(mex_dir)
      candidates <- list.files(root, pattern="tissue_positions.*\\.csv", recursive=TRUE, full.names=TRUE)
      
      # Also check original_root if provided (e.g. tar extraction case)
      if (!is.null(original_root)) {
          candidates <- c(candidates, 
               list.files(original_root, pattern="tissue_positions.*\\.csv", recursive=TRUE, full.names=TRUE)
          )
      }
      
      if (length(candidates) > 0) coords_path <- candidates[1]
  }
  
  coords_df <- NULL
  if (!is.null(coords_path) && file.exists(coords_path)) {
      # Same logic as H5 coord reading
      raw_co <- data.table::fread(coords_path, header=FALSE)
       if (raw_co[1,1] == "barcode") {
          colnames(raw_co) <- as.character(raw_co[1,])
          raw_co <- raw_co[-1,]
      } else {
          colnames(raw_co) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
      }
      raw_co <- raw_co[raw_co$in_tissue == 1, ]
      coords_df <- data.frame(
          barcode = raw_co$barcode,
          imagerow = as.numeric(raw_co$pxl_row),
          imagecol = as.numeric(raw_co$pxl_col)
      )
      
      common <- intersect(colnames(counts), coords_df$barcode)
      counts <- counts[, common, drop=FALSE]
      coords_df <- coords_df[match(common, coords_df$barcode), ]
  }
  
  return(list(counts=counts, coords=coords_df))
}
