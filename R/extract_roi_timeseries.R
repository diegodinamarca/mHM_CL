#' Extract ROI mean time series
#'
#' Computes the spatial mean of a variable inside a region of interest (ROI)
#' for every time step of `mHM_Fluxes_States.nc`.
#'
#' @param nc_path Path to `mHM_Fluxes_States.nc`.
#' @param var_name Name of the variable to extract.
#' @param roi_file Path to a vector file representing the ROI.
#' @param out_file Optional path to write the resulting CSV. If `NULL`, the
#'   data frame is returned invisibly.
#'
#' @return A data frame with columns `time` and `mean`.
#' @examples
#' extract_roi_timeseries("domain/OUT/mHM_Fluxes_States.nc",
#'                        var_name = "Q",
#'                        roi_file = "domain/roi.shp")
#' @export
extract_roi_timeseries <- function(nc_path, var_name, roi_file, out_file = NULL) {
  library(terra)

  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }
  if (!file.exists(roi_file)) {
    stop("ROI file not found: ", roi_file)
  }

  r <- rast(nc_path, subds = var_name)
  roi <- vect(roi_file)
  r <- mask(r, roi)

  means <- as.numeric(global(r, fun = "mean", na.rm = TRUE))
  tvec <- terra::time(r)
  df <- data.frame(time = tvec, mean = means)

  if (!is.null(out_file)) {
    write.csv(df, out_file, row.names = FALSE)
  }

  invisible(df)
}
