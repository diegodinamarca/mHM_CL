#' Get bounding box polygon of an ROI
#'
#' Constructs a polygon from the bounding box of the provided
#' region of interest (ROI).
#'
#' @param roi An [`sf`] object representing the ROI.
#'
#' @return An [`sf`] polygon representing the ROI extent in EPSG:4326.
#' @keywords internal
get_extent <- function(roi) {
  box = st_bbox(roi)
  box.m = as.matrix(rbind(
    c(box[[1]], box[[2]]),
    c(box[[1]], box[[4]]),
    c(box[[3]], box[[4]]),
    c(box[[3]], box[[2]]),
    c(box[[1]], box[[2]])
  ))
  box.pol = st_polygon(list(box.m)) %>% st_sfc(crs = 4326)
}

#' Iterative focal mean filter
#'
#' Applies a 3x3 mean filter recursively `n` times.
#'
#' @param r A [`terra::rast`] object to smooth.
#' @param n Integer. Number of iterations to apply.
#'
#' @return Smoothed [`terra::rast`] object.
#' @keywords internal
focal_repeat = function(r, n){
  if (n != 0){
    r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
    return(focal_repeat(r, n-1))
  }else{
    return(r)
  }
}

#' Extract header from an ESRI ASCII file
#'
#' Reads the first six lines of a `.asc` raster file and writes them
#' to `header_morph.txt` in the specified output folder.
#'
#' @param asc_file Path to the input `.asc` file.
#' @param output_folder Directory where the header file will be created.
#'
#' @return Invisibly returns `NULL`. The header file is written to disk.
#' @keywords internal
extract_asc_header <- function(asc_file, output_folder) {
  # Check if input file exists
  if (!file.exists(asc_file)) {
    stop("The input .asc file does not exist.")
  }
  
  # Check if output folder exists, if not, create it
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Read the first 6 lines of the .asc file (header)
  header_lines <- readLines(asc_file, n = 6)
  
  # Define the path to the output header file
  header_file <- file.path(output_folder, "header_morph.txt")
  
  # Write the header to the file
  writeLines(header_lines, con = header_file)
  
  cat("Header successfully written to", header_file, "\n")
}
