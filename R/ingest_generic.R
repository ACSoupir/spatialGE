#' Ingest Generic Data (DataFrames, Lists, Text Files)
#' @param source InputSource object of type 'generic' or 'list'
#' @export
dispatch_ingest.source_generic <- function(source) {
  # Text files (csv/tsv)
  rna_path <- source$rna
  coord_path <- source$coords
  
  # Read Counts
  counts_df <- data.table::fread(rna_path, header=TRUE, data.table=FALSE)
  # Detect generic vs CosMx logic
  # CosMx usually has "FOV" and "CellID" columns?
  # If so, might need special handling.
  
  # Heuristic: Check if structure is (Gene, Sample1, Sample2...) or (CellID, Gene1, Gene2...) or CosMx wide format?
  # STList generic assumes: Gene in col1 (or rownames), cells in columns.
  # UNLESS it's CosMx text file which is handled differently in legacy STlist?
  # Let's check if it's CosMx-like
  is_cosmx <- any(grepl("fov", colnames(counts_df), ignore.case=TRUE)) && any(grepl("cell_id", colnames(counts_df), ignore.case=TRUE))
  
  if (is_cosmx) {
      # CosMx ingestion logic (splitting by FOV)
      # Adapted from legacy import_smi
      
      # Determine sample/slide name
      slidename <- source$samples
      if (is.null(slidename)) slidename <- basename(rna_path) # Fallback
      
      # Read Coords first (to match FOVs) if present
      spotcoords_df <- NULL
      if (!is.null(coord_path) && file.exists(coord_path)) {
          spotcoords_df <- data.table::fread(coord_path, data.table=FALSE)
      }
      
      results_batch <- list()
      fovid_list <- unique(counts_df[['fov']])
      
      for(fovid in fovid_list) {
          # Subset Counts (Cell x Gene)
          # Legacy filters cell_ID != 0
          sub_df <- counts_df[counts_df[['fov']] == fovid & counts_df[['cell_ID']] != 0, ]
          
          # Construct libnames: fov_X_cell_Y
          libnames <- paste0('fov_', sub_df[['fov']], '_', 'cell_', sub_df[['cell_ID']])
          
          # Remove metadata cols to get count matrix
          # Legacy ignores: fov, cell_ID, cell
          count_mat_df <- sub_df[, !grepl('^fov$|^cell_ID$|^cell$', colnames(sub_df)), drop=FALSE]
          
          # Transpose to Gene x Cell
          # counts_df is Cell x Gene (Genes are columns)
          # So we transpose
          count_mat <- t(as.matrix(count_mat_df))
          colnames(count_mat) <- libnames
          
          # Sparse Matrix
          sp_counts <- Matrix::Matrix(count_mat, sparse=TRUE)
          
          # Coords
          sp_coords <- NULL
          if (!is.null(spotcoords_df)) {
              sub_co <- spotcoords_df[spotcoords_df[['fov']] == fovid, ]
              
              # Construct libnames
              libnames_co <- paste0('fov_', sub_co[['fov']], '_', 'cell_', sub_co[['cell_ID']])
              
              # Select X/Y
              # Legacy expects CenterX_local_px, CenterY_local_px
              # Map to xpos, ypos (Legacy order: libname, xpos, ypos)
              sp_coords <- data.frame(
                  barcode = libnames_co,
                  xpos = as.numeric(sub_co[['CenterX_local_px']]),
                  ypos = as.numeric(sub_co[['CenterY_local_px']])
              )
              
              # Align
              common <- intersect(colnames(sp_counts), sp_coords$barcode)
              sp_counts <- sp_counts[, common, drop=FALSE]
              sp_coords <- sp_coords[match(common, sp_coords$barcode), ]
          }
          
          # Result Name: slidename_fov_X
          res_name <- paste0(slidename, '_fov_', fovid)
          results_batch[[res_name]] <- list(counts=sp_counts, coords=sp_coords)
      }
      
      return(results_batch)
  }
  
  # Assume Col 1 is Genes
  genes <- counts_df[[1]]
  counts <- as.matrix(counts_df[,-1])
  rownames(counts) <- make.unique(as.character(genes))
  counts <- Matrix::Matrix(counts, sparse=TRUE)
  
  # Read Coords
  coords_df <- NULL
  if (!is.null(coord_path) && file.exists(coord_path)) {
      raw_co <- data.table::fread(coord_path, header=TRUE, data.table=FALSE)
      # Expect: Col1=Barcode/CellID, then x/y
      # STList generic: detected col names or assume 1,2,3?
      # Let's look for standard names x/y
      
      # Minimal: barcode, ypos, xpos
      # Just map defaults
      coords_df <- data.frame(
         barcode = as.character(raw_co[[1]]),
         ypos = as.numeric(raw_co[[2]]),
         xpos = as.numeric(raw_co[[3]])
      )
      
      common <- intersect(colnames(counts), coords_df$barcode)
      counts <- counts[, common, drop=FALSE]
      coords_df <- coords_df[match(common, coords_df$barcode), ]
  }
  
  return(list(counts=counts, coords=coords_df))
}

#' @export
dispatch_ingest.source_list <- function(source) {
   # source$rna is a list of data frames (counts)
   # source$coords is a list of data frames (coords) (optional)
   
   # Iterate and Ingest each?
   # But detect_input_source returned ONE source object for the whole list?
   # Actually, yes. `detect_input_source` returns `list(new__source("list", ...))`
   # So this function receives the WHOLE structure.
   
   # STList expects to process this internally.
   # But our Strategy Pattern usually returns (counts=Matrix, coords=DF) for ONE sample.
   # Special Case: 'list' source actually contains MULTIPLE samples.
   # We should probably return a LIST of (counts, coords) objects?
   # Or better: `STList_new` handles the loop over the RESULT of dispatch_ingest?
   
   # If we return a list of results here, `STList_new` needs to know it got a batch.
   # Let's return a NAMED LIST of results.
   
   res <- list()
   rna_list <- source$rna
   coord_list <- source$coords
   
   samples <- names(rna_list)
   for (nm in samples) {
       # Counts
       cnt_df <- rna_list[[nm]]
       # Assume Gene in Col 1
       genes <- cnt_df[[1]]
       vm <- as.matrix(cnt_df[,-1])
       rownames(vm) <- make.unique(as.character(genes))
       sp_counts <- Matrix::Matrix(vm, sparse=TRUE)
       
       # Coords
       sp_coords <- NULL
       if (!is.null(coord_list) && nm %in% names(coord_list)) {
           co_df <- coord_list[[nm]]
           # Assume Barcode, Y, X
           sp_coords <- data.frame(
               barcode = as.character(co_df[[1]]),
               ypos = as.numeric(co_df[[2]]),
               xpos = as.numeric(co_df[[3]])
           )
           
           common <- intersect(colnames(sp_counts), sp_coords$barcode)
           sp_counts <- sp_counts[, common, drop=FALSE]
           sp_coords <- sp_coords[match(common, sp_coords$barcode), ]
       }
       
       res[[nm]] <- list(counts=sp_counts, coords=sp_coords)
   }
   
   return(res)
}

#' Ingest Seurat Object
#' @param source InputSource object of type 'seurat'
#' @export
dispatch_ingest.source_seurat <- function(source) {
  obj <- source$rna
  if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat required")
  
  # Seurat Ingestion (Multi-Sample/Batch Support)
  if (!requireNamespace("Seurat", quietly=TRUE)) stop("Seurat required")
  
  # Helper for GetAssayData (handle V3 vs V5 depreciation)
  # Seurat 5 uses layer='counts', V4/3 uses slot='counts'
  # Try catch or just straightforward call (Seurat warns but works)
  counts_all <- tryCatch({
       Seurat::GetAssayData(obj, slot="counts")
  }, error=function(e) {
       Seurat::GetAssayData(obj, layer="counts")
  })
  
  images <- names(obj@images)
  results_batch <- list()
  
  if (length(images) == 0) {
      # No images? fallback to single object ingestion?
      # STList requires spatial info.
      warning("No images found in Seurat object. Trying global ingestion.")
      # Try entire counts + no coords? Or generic coords?
      # STList typically requires coords.
      # Return single empty coord list?
      return(list(counts=counts_all, coords=NULL))
  }
  
  for (img in images) {
      # 1. Identify Cells in this Slice
      # Seurat::Cells(obj[[img]]) robustly gets global cell IDs for the image
      slice_cells <- tryCatch({
          Seurat::Cells(obj[[img]])
      }, error=function(e) NULL)
      
      if (is.null(slice_cells)) next
      
      # 2. Subset Counts
      common_cells <- intersect(colnames(counts_all), slice_cells)
      if (length(common_cells) == 0) next
      
      sp_counts <- counts_all[, common_cells, drop=FALSE]
      
      # 3. Get Coords
      coords_df <- NULL
      try({
          # Use object-level GetTissueCoordinates for robust column names/IDs
          raw_co <- Seurat::GetTissueCoordinates(obj, image=img)
          
          # Barcodes
          if ("cell" %in% colnames(raw_co)) {
              barcodes <- raw_co$cell
          } else {
              barcodes <- rownames(raw_co)
          }
          
          # Coords X/Y mapping (handle 'imagerow/col' vs 'x/y')
          cn <- colnames(raw_co)
          
          # Logic: Look for x/y or imagecol/imagerow
          # Priority: imagecol/imagerow (Visium standard) > x/y (Generic/Nanostring)
          if ("imagecol" %in% cn && "imagerow" %in% cn) {
              xvec <- raw_co$imagecol
              yvec <- raw_co$imagerow
          } else {
              # Find x/y
              xcol <- cn[grep("^x", cn, ignore.case=TRUE)][1]
              ycol <- cn[grep("^y", cn, ignore.case=TRUE)][1]
              xvec <- raw_co[[xcol]]
              yvec <- raw_co[[ycol]]
          }
          
          coords_df <- data.frame(
              barcode = as.character(barcodes),
              xpos = as.numeric(xvec),
              ypos = as.numeric(yvec)
          )
      }, silent=TRUE)
      
      # Align
      if (!is.null(coords_df)) {
          common <- intersect(colnames(sp_counts), coords_df$barcode)
          sp_counts <- sp_counts[, common, drop=FALSE]
          coords_df <- coords_df[match(common, coords_df$barcode), ]
      }
      
      results_batch[[img]] <- list(counts=sp_counts, coords=coords_df)
  }
  
  return(results_batch)
}
