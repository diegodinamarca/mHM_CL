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

  means <- global(r, fun = "mean", na.rm = TRUE)$mean
  tvec <- terra::time(r)
  df <- data.frame(time = tvec, mean = means)

  if (!is.null(out_file)) {
    write.csv(df, out_file, row.names = FALSE)
  }

  invisible(df)
}

#' Extract ROI mean time series for all variables
#'
#' Computes the spatial mean inside a region of interest (ROI) for
#' every variable available in `mHM_Fluxes_States.nc`.
#'
#' @param nc_path Path to `mHM_Fluxes_States.nc`.
#' @param roi_file Path to a vector file representing the ROI.
#' @param out_file Optional path to write the resulting CSV. If `NULL`, the
#'   data frame is returned invisibly.
#'
#' @return A data frame with a `time` column and one column per variable.
#' @examples
#' extract_roi_timeseries_all("domain/OUT/mHM_Fluxes_States.nc",
#'                            roi_file = "domain/roi.shp")
#' @export
extract_roi_timeseries_all <- function(nc_path, roi_file, out_file = NULL) {
  library(terra)
  library(ncdf4)
  
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }
  if (!file.exists(roi_file)) {
    stop("ROI file not found: ", roi_file)
  }

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  # vars <- c
  vars <- names(nc$var)

  roi <- vect(roi_file)
  all_means <- list()
  time_vec <- NULL
  r <- rast(nc_path)
  time_vec <- as.Date(terra::time(r))
  
  for (v in vars) {
    r <- rast(nc_path, subds = v)
    r <- mask(r, roi)

    means <- global(r, fun = "mean", na.rm = TRUE)$mean
    all_means[[v]] <- means
  }

  df <- data.frame(time = time_vec, all_means)

  if (!is.null(out_file)) {
    write.csv(df, out_file, row.names = FALSE)
  }

  invisible(df)
}

