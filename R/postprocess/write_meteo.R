#' Export aggregated meteorological forcings
#'
#' Reads the daily NetCDF file located in the `meteo` folder of a domain and
#' writes the selected variable either as monthly or yearly layers.
#' Results are saved inside `out/meteo/<var_name>` in the chosen format.
#'
#' @param domain_path Path to the model domain containing the `meteo` folder.
#' @param var_name Name of the variable (and NetCDF file) to export.
#' @param ts Time step of the output. Either "month" or "year".
#' @param roi_mask Logical. Apply the ROI mask to the output.
#' @param roi_file Optional path to a ROI file. Used when `roi_mask` is `TRUE`.
#'   If `NULL`, the `roi_file` defined in `preprocess_config.json` is used.
#' @param out.format Output format: "tif" for GeoTIFF or "nc" for NetCDF.
#'
#' @return Invisibly returns the path of the written file.
#' @examples
#' write_meteo("/path/to/domain", "pre", "year")
#' @export
write_meteo <- function(domain_path, var_name, ts,
                        roi_mask = TRUE, roi_file = NULL,
                        out.format = c("tif", "nc")) {
  library(terra)
  library(jsonlite)

  out.format <- match.arg(out.format)
  ts <- match.arg(tolower(ts), c("month", "year"))

  nc_path <- file.path(domain_path, "meteo", paste0(var_name, ".nc"))
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

  message("Reading variable '", var_name, "' from ", nc_path)
  r <- rast(nc_path, subds = var_name)

  if (roi_mask) {
    if (is.null(roi_file)) {
      config_path <- file.path(domain_path, "preprocess_config.json")
      if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
      }
      config <- fromJSON(config_path)
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    message("Applying ROI mask")
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }

  fun <- if (var_name %in% c("pre", "pet")) "sum" else "mean"
  message("Calculating ", ts, " values")
  if (ts == "month") {
    r <- daily_to_monthly(r, fun = fun)
  } else if (ts == "year") {
    r <- daily_to_monthly(r, fun = fun)
    r <- monthly_to_yearly(r, fun = fun)
  }

  out_dir <- file.path(domain_path, "out", "meteo", var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (out.format == "tif") {
    fname <- file.path(out_dir, paste0(var_name, "_", ts, ".tif"))
    message("Writing output to ", fname)
    writeRaster(r, fname, overwrite = TRUE)
  } else {
    fname <- file.path(out_dir, paste0(var_name, "_", ts, ".nc"))
    message("Writing output to ", fname)
    writeCDF(r, filename = fname, varname = var_name, overwrite = TRUE)
  }

  invisible(fname)
}

#' Export long term annual mean of a meteorological variable
#'
#' Computes the long term annual mean of a daily forcing located in the
#' `meteo` folder. The result is written inside `out/meteo/<var_name>`.
#'
#' @inheritParams write_meteo
#' @return Invisibly returns the path of the written file.
#' @examples
#' write_meteo_annual_mean("/path/to/domain", "pre")
#' @export
write_meteo_annual_mean <- function(domain_path, var_name,
                                    roi_mask = TRUE, roi_file = NULL,
                                    out.format = c("tif", "nc")) {
  library(terra)
  library(jsonlite)

  out.format <- match.arg(out.format)

  nc_path <- file.path(domain_path, "meteo", paste0(var_name, ".nc"))
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

  message("Reading variable '", var_name, "' from ", nc_path)
  r <- rast(nc_path, subds = var_name)

  if (roi_mask) {
    if (is.null(roi_file)) {
      config_path <- file.path(domain_path, "preprocess_config.json")
      if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
      }
      config <- fromJSON(config_path)
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    message("Applying ROI mask")
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }

  fun <- if (var_name %in% c("pre", "pet")) "sum" else "mean"
  message("Calculating annual mean")
  r_month <- daily_to_monthly(r, fun = fun)
  r_ann <- annual_mean(r_month, fun = fun)

  out_dir <- file.path(domain_path, "out", "meteo", var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (out.format == "tif") {
    fname <- file.path(out_dir, paste0(var_name, "_annual.tif"))
    message("Writing output to ", fname)
    writeRaster(r_ann, fname, overwrite = TRUE)
  } else {
    fname <- file.path(out_dir, paste0(var_name, "_annual.nc"))
    message("Writing output to ", fname)
    writeCDF(r_ann, filename = fname, varname = var_name, overwrite = TRUE)
  }

  invisible(fname)
}

#' Process daily meteorological data to multiple time scales
#'
#' Convenience workflow that exports monthly, yearly and annual mean rasters
#' for a meteorological variable.
#'
#' @inheritParams write_meteo
#' @return Invisibly returns the written file paths as a named vector.
#' @examples
#' process_meteo_variable("/path/to/domain", "pre")
#' @export
process_meteo_variable <- function(domain_path, var_name,
                                   roi_mask = TRUE, roi_file = NULL,
                                   out.format = c("tif", "nc")) {
  library(terra)
  library(jsonlite)

  out.format <- match.arg(out.format)

  nc_path <- file.path(domain_path, "meteo", paste0(var_name, ".nc"))
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

  message("Reading variable '", var_name, "' from ", nc_path)
  r <- rast(nc_path, subds = var_name)

  if (roi_mask) {
    if (is.null(roi_file)) {
      config_path <- file.path(domain_path, "preprocess_config.json")
      if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
      }
      config <- fromJSON(config_path)
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    message("Applying ROI mask")
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }

  fun <- if (var_name %in% c("pre", "pet")) "sum" else "mean"
  message("Calculating monthly values")
  r_month <- daily_to_monthly(r, fun = fun)
  message("Calculating yearly values")
  r_year  <- monthly_to_yearly(r_month, fun = fun)
  message("Calculating long term annual mean")
  r_ann   <- annual_mean(r_month, fun = fun)

  out_dir <- file.path(domain_path, "out", "meteo", var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (out.format == "tif") {
    fname_month <- file.path(out_dir, paste0(var_name, "_month.tif"))
    fname_year  <- file.path(out_dir, paste0(var_name, "_year.tif"))
    fname_ann   <- file.path(out_dir, paste0(var_name, "_annual.tif"))
    message("Writing monthly raster to ", fname_month)
    writeRaster(r_month, fname_month, overwrite = TRUE)
    message("Writing yearly raster to ", fname_year)
    writeRaster(r_year,  fname_year,  overwrite = TRUE)
    message("Writing annual mean raster to ", fname_ann)
    writeRaster(r_ann,   fname_ann,   overwrite = TRUE)
  } else {
    fname_month <- file.path(out_dir, paste0(var_name, "_month.nc"))
    fname_year  <- file.path(out_dir, paste0(var_name, "_year.nc"))
    fname_ann   <- file.path(out_dir, paste0(var_name, "_annual.nc"))
    message("Writing monthly raster to ", fname_month)
    writeCDF(r_month, filename = fname_month, varname = var_name, overwrite = TRUE)
    message("Writing yearly raster to ", fname_year)
    writeCDF(r_year,  filename = fname_year,  varname = var_name, overwrite = TRUE)
    message("Writing annual mean raster to ", fname_ann)
    writeCDF(r_ann,   filename = fname_ann,   varname = var_name, overwrite = TRUE)
  }

  invisible(c(monthly = fname_month, yearly = fname_year, annual = fname_ann))
}
