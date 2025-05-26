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

focal_repeat = function(r, n){
  if (n != 0){
    r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
    return(focal_repeat(r, n-1))
  }else{
    return(r)
  }
}

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