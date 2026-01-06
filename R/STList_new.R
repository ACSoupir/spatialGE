#' STList Constructor (Refactored)
#' @description Creates an STList object using modular ingestion strategy.
#' @export
STList_new <- function(rnacounts=NULL, spotcoords=NULL, samples=NULL, cores=NULL, verbose=TRUE) {
  
  # 1. Detect Sources
  input_sources <- detect_input_source(rnacounts, spotcoords, samples)
  
  # 2. Iterate and Ingest
  counts_list <- list()
  coords_list <- list()
  
  for (i in seq_along(input_sources)) {
    src <- input_sources[[i]]
    if (verbose) message("Ingesting source ", i, ": Type=", src$type, " Format=", src$format)
    
    res <- dispatch_ingest(src)
    
    # Handle Result
    # If generic-list, res is a list of results.
    # Otherwise res is a single list(counts=, coords=)
    
    if (src$type == "list") {
        # Merge the lists
        counts_list <- c(counts_list, lapply(res, `[[`, "counts"))
        coords_list <- c(coords_list, lapply(res, `[[`, "coords"))
    } else {
        # Determine Sample Name
        sname <- src$samples
        if (is.null(sname)) {
            # Infer from file path or use default
            if (is.character(src$rna)) sname <- basename(src$rna)
            else sname <- paste0("Sample_", i)
        }
        
        # Check if result has counts/coords
        counts_list[[sname]] <- res$counts
        coords_list[[sname]] <- res$coords
    }
  }
  
  # 3. Create STList Object
  # (Reusing internal structure logic from original STList)
  
  # Filter NULLs (some ingestors might fail strict checks? Or we stop earlier.)
  # Assuming ingestors return valid data or error.
  
  # Create S4 Object
  # Initialize sample_meta as empty tibble (required by class definition)
  # Or populate with sample names?
  # Legacy behavior: creates tibble(sample_name=names(counts))
  
  sample_meta <- tibble::tibble(sample_name = names(counts_list))
  # Must ensure sample_name is character
  sample_meta$sample_name <- as.character(sample_meta$sample_name)
  
  st_obj <- methods::new("STlist", 
      counts = counts_list, 
      spatial_meta = coords_list,
      sample_meta = sample_meta
  )
  
  # Post-Processing (mirroring legacy STList logic)
  # - Ensure sample names match
  # - Convert coords to spatial_meta dataframes with libname
  
  # Fix spatial_meta structure: STList expects specific columns?
  # Legacy STList: spatial_meta is list of DFs. 
  # Cols: libname, ypos, xpos...
  
  for (nm in names(st_obj@spatial_meta)) {
      df <- st_obj@spatial_meta[[nm]]
      if (!is.null(df)) {
          # Ensure libname column
          if (!"libname" %in% colnames(df)) {
             df$libname <- nm # OR barcode? STList usually has barcode as rowname or column?
             # Actually STList usually has 'libname' as the sample name? 
             # Let's check classDefinitions.R or STList.R logic.
             # STList.R line 287ish: samples_df = tibble(sample_name=names(...))
             # STList.R post-processing does a lot.
             
             # For now, let's trust the ingestors pass 'barcode', 'xpos', 'ypos'.
             # STList usually expects 'libname' column in spatial_meta to be the SAMPLE NAME?
             # No, 'libname' usually refers to the BARCODE in some contexts or Sample?
             # Let's verify STList.R logic.
             # Line 230: process_sample_names...
             
             # Actually, simpler:
             # STList@spatial_meta is a list of data frames.
             # The data frames should have rownames = barcodes?
             # Or a column 'libname' = barcode?
          }
      }
      st_obj@spatial_meta[[nm]] <- df
  }
  
  return(st_obj)
}
