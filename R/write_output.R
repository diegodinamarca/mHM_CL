#' Export mHM outputs by time step
#'
#' Reads `mHM_Fluxes_States.nc` located in the `OUT` folder of a domain and
#' writes the selected variable either as monthly or yearly layers. Results are
#' saved inside `OUT/<var_name>` in the chosen format.
#'
#' @param domain_path Path to the model domain containing the `OUT` folder.
#' @param var_name Name of the variable inside `mHM_Fluxes_States.nc` to export.
#' @param ts Time step of the output. Either "month" or "year".
#' @param roi_mask Logical. Apply the ROI mask to the output.
#' @param roi_file Optional path to a ROI file. Used when `roi_mask` is `TRUE`.
#'   If `NULL`, the `roi_file` defined in `preprocess_config.json` is used.
#' @param out.format Output format: "tif" for GeoTIFF or "nc" for NetCDF.
#'
#' @return Invisibly returns the paths of the written files.
#' @examples
#' write_output("/path/to/domain", "Q", "year")
#' @export
write_output <- function(domain_path, var_name, ts,
                         roi_mask = TRUE, roi_file = NULL,
                         out.format = c("tif", "nc")) {
  library(terra)
  library(jsonlite)

  out.format <- match.arg(out.format)
  ts <- match.arg(tolower(ts), c("month", "year"))

  nc_path <- file.path(domain_path, "OUT", "mHM_Fluxes_States.nc")
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

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
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }

  if (ts == "month") {
    r <- daily_to_monthly(r, fun = "mean")
  } else if (ts == "year") {
    r <- monthly_to_yearly(r, fun = "mean")
  }

  out_dir <- file.path(domain_path, "OUT", var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (out.format == "tif") {
    fname <- file.path(out_dir, paste0(var_name, ".tif"))
    writeRaster(r, fname, overwrite = TRUE)
    paths <- fname
    paths <- vector("character", nlyr(r))
    for (i in seq_len(nlyr(r))) {
      suffix <- if (ts == "month") "%Y_%m" else "%Y"
      dstr <- format(time(r)[i], suffix)
      fname <- file.path(out_dir, paste0(var_name, "_", dstr, ".tif"))
      writeRaster(r[[i]], fname, overwrite = TRUE)
      paths[i] <- fname
    }
  } else {
    fname <- file.path(out_dir, paste0(var_name, ".nc"))
    writeCDF(r, filename = fname, varname = var_name, overwrite = TRUE)
    paths <- fname
  }

  invisible(paths)
}

