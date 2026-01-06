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
      # CosMx ingestion logic (adapted from legacy)
      # Usually long format or wide?
      # For now, implementing standard generic assumption: Gene x Cell matrix
      # If specialized CosMx needed, we can augment this later or create `ingest_cosmx`.
      # Assume standard generic for now.
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
  
  # Extract counts
  # Handle V3 vs V5?
  # SeuratObject::GetAssayData(obj, slot="counts")
  counts <- Seurat::GetAssayData(obj, slot="counts")
  
  # Extract Coords
  # obj@images
  coords_df <- NULL
  if (length(obj@images) > 0) {
      # Use first image?
      img <- obj@images[[1]]
      coords <- Seurat::GetTissueCoordinates(img)
      # Seurat coords: y, x (rows, columns?)
      # Varies by version.
      # Usually: c("imagerow", "imagecol")
      
      coords_df <- data.frame(
          barcode = rownames(coords),
          ypos = coords[,1], # verify mapping order?
          xpos = coords[,2]
      )
  }
  
  # If multiple images/samples in Seurat object?
  # STList logic was complex. For now, flattened ingest of whole object.
  # Or slice by image?
  
  return(list(counts=counts, coords=coords_df))
}
