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

#' Mosaic mHM outputs from multiple domains
#'
#' Reads `mHM_Fluxes_States.nc` from each domain, applies the ROI mask defined
#' in the respective `preprocess_config.json` and merges all domains into
#' separate NetCDF files, one per variable.
#'
#' @param domains Character vector with paths to the domain folders. If `NULL`,
#'   all folders in the current working directory matching "domain" are used.
#' @param out_dir Folder where the mosaicked NetCDFs will be created.
#' @param vars Character vector of variable names to mosaic. By default the
#'   variables used in `visualize_annual_outputs()` except precipitation and
#'   potential evapotranspiration are processed.
#' @examples
#' mosaic_outputs(domains = c("domain_1", "domain_2"),
#'                out_dir = "domain_chile/OUT")
#' @export
mosaic_outputs <- function(domains = NULL,
                           out_dir = "domain_chile/OUT",
                           vars = c("snowpack", "SM_Lall", "satSTW", "aET", "Q")) {
  library(terra)
  library(jsonlite)
  
  if (is.null(domains)) {
    domains <- dir(pattern = "domain", full.names = TRUE)
  }
  if (length(domains) == 0) {
    stop("No domain folders found")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  written <- vector("character", length(vars))
  for (v in vars) {
    message("Mosaicking variable: ", v)
    rasters <- list()
    for (d in domains) {
      cfg_path <- file.path(d, "preprocess_config.json")
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- fromJSON(cfg_path)
      nc_file <- file.path(d, cfg$out_folder, "mHM_Fluxes_States.nc")
      if (!file.exists(nc_file)) {
        stop("NetCDF file not found: ", nc_file)
      }
      r <- rast(nc_file, subds = v)
      if (!is.null(cfg$roi_file) && file.exists(cfg$roi_file)) {
        roi <- vect(cfg$roi_file)
        r <- mask(r, roi)
      }
      rasters[[length(rasters) + 1]] <- r
    }
    
    mosaic_r <- do.call(mosaic, rasters)
    out_file <- file.path(out_dir, paste0(v, ".nc"))
    writeCDF(mosaic_r, filename = out_file, varname = v, overwrite = TRUE)
    written[which(vars == v)] <- out_file
  }
  
  invisible(written)
}

#' Mosaic meteorological forcings from multiple domains
#'
#' Reads daily NetCDF files located in the `meteo` folder of each domain and
#' merges them into a single file per variable. The ROI mask defined in each
#' `preprocess_config.json` is applied prior to mosaicking.
#'
#' @param domains Character vector with paths to the domain folders. If `NULL`,
#'   all folders in the current working directory matching "domain" are used.
#' @param out_dir Folder where the mosaicked NetCDFs will be created.
#' @param vars Character vector of forcing names to mosaic. Defaults to
#'   `c("pre", "pet", "tmin", "tmax", "tavg")` which correspond to
#'   precipitation, potential evapotranspiration and temperature metrics.
#' @examples
#' mosaic_meteo(domains = c("domain_1", "domain_2"),
#'              out_dir = "domain_chile/OUT")
#' @export
mosaic_meteo <- function(domains = NULL,
                         out_dir = "domain_chile/OUT",
                         vars = c("pre", "pet", "tmin", "tmax", "tavg")) {
  library(terra)
  library(jsonlite)
  
  if (is.null(domains)) {
    domains <- dir(pattern = "domain", full.names = TRUE)
  }
  if (length(domains) == 0) {
    stop("No domain folders found")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  written <- vector("character", length(vars))
  for (v in vars) {
    message("Mosaicking variable: ", v)
    rasters <- list()
    for (d in domains) {
      cfg_path <- file.path(d, "preprocess_config.json")
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- fromJSON(cfg_path)
      meteo_file <- file.path(d, "meteo", paste0(v, ".nc"))
      if (!file.exists(meteo_file)) {
        stop("NetCDF file not found: ", meteo_file)
      }
      r <- rast(meteo_file, subds = v)
      if (!is.null(cfg$roi_file) && file.exists(cfg$roi_file)) {
        roi <- vect(cfg$roi_file)
        r <- mask(r, roi)
      }
      rasters[[length(rasters) + 1]] <- r
    }
    
    mosaic_r <- do.call(mosaic, rasters)
    out_file <- file.path(out_dir, paste0(v, ".nc"))
    writeCDF(mosaic_r, filename = out_file, varname = v, overwrite = TRUE)
    written[which(vars == v)] <- out_file
  }
  
  invisible(written)
}


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
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".tif"))
    writeRaster(r, fname, overwrite = TRUE)
    paths <- fname
    paths <- vector("character", nlyr(r))
  } else {
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".nc"))
    writeCDF(r, filename = fname, varname = var_name, overwrite = TRUE)
    paths <- fname
  }
  
  invisible(paths)
}

#' Visualize annual mean model outputs and forcings
#'
#' Reads the standard mHM NetCDF outputs together with precipitation and
#' potential evapotranspiration forcings and displays the long term annual means
#' in a 2x4 panel plot.
#'
#' @param domain_path Path to the model domain containing the `OUT` and `meteo`
#'   folders.
#' @param mask_roi Logical. If `TRUE`, the outputs are masked using the
#'   `roi_file` specified in `preprocess_config.json` located inside
#'   `domain_path`.
#' @examples
#' visualize_annual_outputs("/path/to/domain")
#' @export
visualize_annual_outputs <- function(domain_path, mask_roi = FALSE){
  library(terra)
  library(jsonlite)
  roi <- NULL
  if (mask_roi) {
    config_path <- file.path(domain_path, "preprocess_config.json")
    if (!file.exists(config_path)) {
      stop("Configuration file not found: ", config_path)
    }
    config <- jsonlite::fromJSON(config_path)
    roi_path <- config$roi_file
    if (!grepl("^(/|[A-Za-z]:)", roi_path)) {
      roi_path <- file.path(roi_path)
    }
    if (!file.exists(roi_path)) {
      stop("ROI file not found: ", roi_path)
    }
    roi <- vect(roi_path)
  }
  out_file <- file.path(domain_path, "OUT", "mHM_Fluxes_States.nc")
  meteo_folder <- file.path(domain_path, "meteo")
  
  snowpack <- rast(out_file, subds = "snowpack")
  soilmoist <- rast(out_file, subds = "SM_Lall")
  satstw <- rast(out_file, subds = "satSTW")
  aet <- rast(out_file, subds = "aET")
  runoff <- rast(out_file, subds = "Q")
  
  pre <- rast(file.path(meteo_folder, "pre.nc"), subds = "pre")
  pet <- rast(file.path(meteo_folder, "pet.nc"), subds = "pet")
  
  message("Processing precipitation")
  pre_ann <- annual_mean(pre, "sum")
  message("Processing potential evapotranspiration")
  pet_ann <- annual_mean(pet, "sum")
  message("Processing actual evapotranspiration")
  aet_ann <- annual_mean(aet, "sum")
  message("Processing runoff")
  runoff_ann <- annual_mean(runoff, "sum")
  
  message("Processing soil moisture")
  sm_ann <- annual_mean(soilmoist, "mean")
  message("Processing snowpack")
  snow_ann <- annual_mean(snowpack, "mean")
  message("Processing groundwater level")
  sat_ann <- annual_mean(satstw, "mean")
  if (!is.null(roi)) {
    pre_ann <- mask(pre_ann, roi)
    pet_ann <- mask(pet_ann, roi)
    aet_ann <- mask(aet_ann, roi)
    runoff_ann <- mask(runoff_ann, roi)
    disp_ann = pre_ann - aet_ann
    
    sm_ann <- mask(sm_ann, roi)
    snow_ann <- mask(snow_ann, roi)
    sat_ann <- mask(sat_ann, roi)
  }
  
  flux_range <- range(c(values(pre_ann), values(pet_ann),
                        values(aet_ann), values(runoff_ann)), na.rm = TRUE)
  state_range <- range(c(values(snow_ann), values(sm_ann), values(sat_ann)),
                       na.rm = TRUE)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 5), mar = c(4, 4, 2, 5))
  
  plot(pre_ann, main = "Precipitation", zlim = flux_range)
  plot(pet_ann, main = "Potential ET", zlim = flux_range)
  plot(aet_ann, main = "Actual ET", zlim = flux_range)
  plot(runoff_ann, main = "Runoff", zlim = flux_range)
  plot(disp_ann, main = "Pr-ET", zlim = flux_range)
  
  plot(sm_ann, main = "Soil moisture", zlim = state_range)
  plot(snow_ann, main = "Snowpack", zlim = state_range)
  plot(sat_ann, main = "GW level", zlim = state_range)
  plot.new()
}