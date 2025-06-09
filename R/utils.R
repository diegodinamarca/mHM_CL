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

#' Compute annual mean across years
#'
#' Summarizes a SpatRaster with monthly layers by first aggregating to yearly
#' values using either the sum or mean of the months, and then averaging across
#' all available years.
#'
#' @param r A [`terra::SpatRaster`] with a monthly time dimension.
#' @param fun Character. Either "mean" or "sum" specifying how to aggregate the
#'   months of each year.
#' @return A single-layer SpatRaster representing the long term annual mean.
#' @export
annual_mean <- function(r, fun = c("mean", "sum")) {
  fun <- match.arg(fun)
  tvec <- terra::time(r)
  years <- format(tvec, "%Y")
  idx <- split(seq_along(years), years)
  yearly <- lapply(idx, function(i) {
    if (fun == "sum") {
      sum(r[[i]], na.rm = TRUE)
    } else {
      mean(r[[i]], na.rm = TRUE)
    }
  })
  yearly <- rast(yearly)
  mean(yearly, na.rm = TRUE)
}

#' Aggregate daily data to monthly values
#'
#' Converts a SpatRaster with daily layers into monthly sums or means depending
#' on `fun`.
#'
#' @param r A [`terra::SpatRaster`] with a daily time dimension.
#' @param fun Character. Either "mean" or "sum" specifying how to aggregate the
#'   days of each month.
#' @return A SpatRaster with one layer per month.
#' @export
daily_to_monthly <- function(r, fun = c("mean", "sum")) {
  fun <- match.arg(fun)
  tvec <- terra::time(r)
  ym <- format(tvec, "%Y-%m")
  idx <- split(seq_along(ym), ym)
  monthly <- lapply(idx, function(i) {
    if (fun == "sum") {
      sum(r[[i]], na.rm = TRUE)
    } else {
      mean(r[[i]], na.rm = TRUE)
    }
  })
  monthly <- rast(monthly)
  time(monthly) <- as.Date(paste0(names(idx), "-01"))
  monthly
}

#' Aggregate monthly data to yearly values
#'
#' Converts a SpatRaster with monthly layers into yearly sums or means depending
#' on `fun`.
#'
#' @param r A [`terra::SpatRaster`] with a monthly time dimension.
#' @param fun Character. Either "mean" or "sum" specifying how to aggregate the
#'   months of each year.
#' @return A SpatRaster with one layer per year.
#' @export
monthly_to_yearly <- function(r, fun = c("mean", "sum")) {
  fun <- match.arg(fun)
  tvec <- terra::time(r)
  years <- format(tvec, "%Y")
  idx <- split(seq_along(years), years)
  yearly <- lapply(idx, function(i) {
    if (fun == "sum") {
      sum(r[[i]], na.rm = TRUE)
    } else {
      mean(r[[i]], na.rm = TRUE)
    }
  })
  yearly <- rast(yearly)
  time(yearly) <- as.Date(paste0(names(idx), "-01-01"))
  yearly
}
