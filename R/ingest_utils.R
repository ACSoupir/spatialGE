#' @title InputSource Class Definitions
#' @description Defines structure for various input sources detected by DetectInput
#' @keywords internal
#' @name ingest_utils
NULL

#' Constructor for InputSource
#' @param type character string of input type (visium, xenium, generic, seurat, list)
#' @param format character string detailed format (h5, mex, csv, etc)
#' @param rna object or path for counts
#' @param coords object or path for coordinates
#' @param samples object or path for samples metadata
#' @return list with class attribute set
new_input_source <- function(type, format=NULL, rna=NULL, coords=NULL, samples=NULL) {
  structure(
    list(
      type = type,
      format = format,
      rna = rna,
      coords = coords,
      samples = samples
    ),
    class = c(paste0("source_", type), "InputSource")
  )
}

#' @title Detect Input Source
#' @description Analyzes arguments to determine the Input Source strategy
detect_input_source <- function(rnacounts=NULL, spotcoords=NULL, samples=NULL) {
  
  # 1. Seurat Object
  if (inherits(rnacounts, "Seurat")) {
    return(list(new_input_source("seurat", rna=rnacounts, samples=samples)))
  }
  
  # 2. List of DataFrames
  if (inherits(rnacounts, "list") && !is.null(names(rnacounts)) && 
      all(sapply(rnacounts, is.data.frame))) {
    return(list(new_input_source("list", rna=rnacounts, coords=spotcoords, samples=samples)))
  }
  
  # 3. File Paths (Vector)
  # Pre-check: if rnacounts is a char vector, we iterate or treat as batch
  if (is.character(rnacounts)) {
    sources <- list()
    for (i in seq_along(rnacounts)) {
      path <- rnacounts[i]
      
      # Associated sample name
      current_sample <- if (!is.null(samples) && length(samples) == length(rnacounts)) samples[i] else NULL
      
      # Associated coord file
      current_coord <- if (!is.null(spotcoords) && length(spotcoords) == length(rnacounts)) spotcoords[i] else NULL
      
      # Check Directory (Visium / Xenium / GeoMx)
      if (dir.exists(path)) {
        # Check for Xenium specific file (cell_feature_matrix)
        # Check for Visium specific file (filtered_feature_bc_matrix)
        
        # Helper to look for file patterns
        has_file <- function(p, pattern) length(list.files(p, pattern=pattern, recursive=TRUE)) > 0
        
        if (has_file(path, "cell_feature_matrix")) {
           sources[[i]] <- new_input_source("xenium", rna=path, coords=current_coord, samples=current_sample)
        } else if (has_file(path, "filtered_feature_bc_matrix") || has_file(path, "raw_feature_bc_matrix") ||
                   has_file(path, "tissue_positions")) {
           sources[[i]] <- new_input_source("visium", rna=path, coords=current_coord, samples=current_sample)
        } else if (has_file(path, ".dcc$")) {
           sources[[i]] <- new_input_source("geomx", rna=path, samples=samples) # GeoMx usually takes directory
        } else {
           # Fallback or error?
           # Maybe treat as generic directory but warn?
           warning("Directory ", path, " recognized but no known spatial format found inside. Treating as generic.")
           sources[[i]] <- new_input_source("generic", rna=path, coords=current_coord, samples=current_sample)
        }
      } 
      
      # Check File (H5 / CSV / TSV)
      else if (file.exists(path)) {
        if (grepl("\\.h5$", path)) {
            # H5 could be Visium or Xenium or Generic
            # Usually strict separation requires checking content, but simpler heuristic:
            # If spotcoords provided, likely generic H5 or specific overrides
            # If not provided, usually 10x H5 containers have both
            # Let's verify h5 file structure inside ingestor, assume '10x_h5' for now if generic
            sources[[i]] <- new_input_source("h5_10x", rna=path, coords=current_coord, samples=current_sample)
        } else {
            # Text file (CSV/TSV/CosMx)
            sources[[i]] <- new_input_source("generic", rna=path, coords=current_coord, samples=current_sample)
        }
      }
    }
    return(sources)
  }
  
  stop("Unknown input format")
}

#' @title Dispatch Ingest
#' @description Generic function to dispatch parsing logic based on InputSource class
dispatch_ingest <- function(source) {
  UseMethod("dispatch_ingest", source)
}

#' Default dispatch logic
dispatch_ingest.default <- function(source) {
  stop("No ingestor defined for input type: ", source$type)
}

#' Safe Matrix Reader
#' @description Robustly reads Matrix market files, handling .gz via temp expansion
#' @keywords internal
read_mtx_safe <- function(path) {
  if (grepl("\\.gz$", path)) {
      # Decompress to temporary file
      tmp_out <- tempfile(fileext = ".mtx")
      
      con_in <- gzfile(path, open="rb")
      con_out <- file(tmp_out, open="wb")
      on.exit({
         try(if (isValid(con_in)) close(con_in), silent=TRUE)
         # isValid? No such function base. try(close) is enough.
         try(close(con_in), silent=TRUE)
         try(close(con_out), silent=TRUE)
         if (file.exists(tmp_out)) unlink(tmp_out)
      })
      
      # Copy stream
      while(length(buf <- readBin(con_in, raw(), n=1024*1024))) {
          writeBin(buf, con_out)
      }
      
      # Flush/Close outputs explicitly to ensure write completion before read
      try(close(con_out), silent=TRUE)
      try(close(con_in), silent=TRUE) # explicit close for tidiness
      
      # Read text matrix from TEMP file (uncompressed)
      return(Matrix::readMM(tmp_out))
  } else {
      return(Matrix::readMM(path))
  }
}
