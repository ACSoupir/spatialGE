#' Ingest Xenium Data
#' @param source InputSource object of type 'xenium'
#' @export
dispatch_ingest.source_xenium <- function(source) {
  path <- source$rna
  
  # 1. Look for H5
  h5_files <- list.files(path, pattern = "cell_feature_matrix.h5$", full.names = TRUE, recursive = TRUE)
  if (length(h5_files) > 0) {
      return(ingest_xenium_h5(h5_files[1], source$coords, source$samples))
  }
  
  # 2. Look for MEX (.tar.gz or folder)
  mex_dirs <- list.files(path, pattern = "cell_feature_matrix$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  mex_tars <- list.files(path, pattern = "cell_feature_matrix.tar.gz$", full.names = TRUE, recursive = TRUE)
  
  if (length(mex_dirs) > 0) {
      return(ingest_xenium_mex_dir(mex_dirs[1], source$coords, source$samples))
  } else if (length(mex_tars) > 0) {
      tmp_ex <- file.path(tempdir(), paste0("xen_ex_", basename(path)))
      utils::untar(mex_tars[1], exdir=tmp_ex)
      inner_dirs <- list.files(tmp_ex, full.names=TRUE, include.dirs=TRUE)
      # Pass 'path' (original source root) as backup for coords lookup
      return(ingest_xenium_mex_dir(inner_dirs[1], source$coords, source$samples, original_root=path))
  }
  
  stop("Could not find valid Xenium H5 or MEX data in: ", path)
}

#' Internal Xenium H5 Reader
ingest_xenium_h5 <- function(h5_path, coords_path=NULL, sample_name=NULL) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("hdf5r required")
  
  h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(h5_file$close_all())
  
  root <- names(h5_file)[1] # usually 'matrix'
  grp <- h5_file[[root]]
  
  data <- grp[["data"]][]
  indices <- grp[["indices"]][] 
  indptr <- grp[["indptr"]][]
  shape <- grp[["shape"]][]
  
  # Xenium features usually strictly named 'features' in recent outputs
  feat_grp <- grp[["features"]]
  genes <- feat_grp[["name"]][]
  barcodes <- grp[["barcodes"]][]
  
  counts <- Matrix::sparseMatrix(
      x = as.numeric(data),
      i = indices,
      p = indptr,
      dims = shape,
      index1 = FALSE
  )
  rownames(counts) <- make.unique(genes)
  colnames(counts) <- barcodes
  
  # Coordinates (cells.parquet or cells.csv)
  if (is.null(coords_path)) {
      # Try known patterns in parent dir
      parent <- dirname(dirname(h5_path))
      candidates <- c(
          file.path(parent, "cells.parquet"),
          file.path(parent, "cells.csv")
      )
      valid <- candidates[file.exists(candidates)]
      if (length(valid) > 0) coords_path <- valid[1]
  }
  
  coords_df <- NULL
  if (!is.null(coords_path) && file.exists(coords_path)) {
      if (grepl("\\.parquet$", coords_path)) {
          if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow required for parquet")
          raw_co <- arrow::read_parquet(coords_path)
      } else {
          raw_co <- data.table::fread(coords_path)
      }
      
      # Xenium standard: cell_id, x_centroid, y_centroid
      # Map to STList: barcode, xpos, ypos (Warning: STList uses 'imagerow/col', but xenium is usually physical micron coords?)
      # STList usually expects 'imagerow'/'imagecol'. 
      # Let's check original import_xenium.R logic...
      # Original mapped: barcode=cell_id, ypos=y_centroid, xpos=x_centroid
      # Note: STlist structure needs 'ypos' and 'xpos' in @spatial_meta? Or coords DF?
      # Original import_xenium returned list(rawcounts=..., coords=DF(barcode, ypos, xpos))
      # Let's stick to that.
      
      coords_df <- data.frame(
          barcode = raw_co$cell_id,
          ypos = raw_co$y_centroid,
          xpos = raw_co$x_centroid
      )
      
      # Optional: Add imagerow/col aliases if STList strictness requires it?
      # STList::process_sample_names -> ... -> STList checks 'imagerow' presence?
      # Actually, detecting spatial Platform in STList depends on columns.
      # But standardizing to 'ypos'/'xpos' is cleaner for Xenium.
      
      common <- intersect(colnames(counts), coords_df$barcode)
      counts <- counts[, common, drop=FALSE]
      coords_df <- coords_df[match(common, coords_df$barcode), ]
  }
  
  return(list(counts=counts, coords=coords_df))
}

#' Internal Xenium MEX Directory Reader
ingest_xenium_mex_dir <- function(mex_dir, coords_path=NULL, sample_name=NULL, original_root=NULL) {
  # Logic: Read matrix.mtx, features.tsv, barcodes.tsv
  # Same as Visium MEX but possibly different feature file structure?
  # Usually standard MTX.
  
  mtx_path <- file.path(mex_dir, "matrix.mtx")
  if (!file.exists(mtx_path)) {
      mtx_path <- paste0(mtx_path, ".gz")
  }
  
  if (!file.exists(mtx_path)) stop("Matrix file not found: ", mex_dir)
  
  counts <- read_mtx_safe(mtx_path)
  
  feat_path <- file.path(mex_dir, "features.tsv")
  if (!file.exists(feat_path)) feat_path <- paste0(feat_path, ".gz")
  genes <- data.table::fread(feat_path, header=FALSE)[[2]]
  rownames(counts) <- make.unique(genes)
  
  bc_path <- file.path(mex_dir, "barcodes.tsv")
  if (!file.exists(bc_path)) bc_path <- paste0(bc_path, ".gz")
  bcs <- data.table::fread(bc_path, header=FALSE)[[1]]
  colnames(counts) <- bcs
  
  # Coordinates
  if (is.null(coords_path)) {
      # Try candidates in mex_dir parent (default) OR original_root
      candidates <- c()
      parent <- dirname(mex_dir)
      
      # 1. Parent of mex dir (if unzipped in place)
      candidates <- c(candidates, 
          file.path(parent, "cells.parquet"),
          file.path(parent, "cells.csv")
      )
      
      # 2. Original root (if tarball extracted to temp)
      if (!is.null(original_root)) {
          candidates <- c(candidates,
              file.path(original_root, "cells.parquet"),
              file.path(original_root, "cells.csv")
          )
      }
      
      valid <- candidates[file.exists(candidates)]
      if (length(valid) > 0) coords_path <- valid[1]
  }
  
  coords_df <- NULL
  if (!is.null(coords_path) && file.exists(coords_path)) {
       if (grepl("\\.parquet$", coords_path)) {
          if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow required for parquet")
          raw_co <- arrow::read_parquet(coords_path)
      } else {
          raw_co <- data.table::fread(coords_path)
      }
      coords_df <- data.frame(
          barcode = raw_co$cell_id,
          ypos = raw_co$y_centroid,
          xpos = raw_co$x_centroid
      )
      common <- intersect(colnames(counts), coords_df$barcode)
      counts <- counts[, common, drop=FALSE]
      coords_df <- coords_df[match(common, coords_df$barcode), ]
  }
  
  return(list(counts=counts, coords=coords_df))
}
