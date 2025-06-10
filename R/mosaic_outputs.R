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

#' Visualize annual mean meteorological forcings across domains
#'
#' Instead of writing full mosaicked NetCDFs, this function computes the long
#' term annual mean of each forcing and displays them in a multipanel plot.
#'
#' @param domains Character vector with paths to domain folders. If `NULL`, all
#'   folders matching "domain" in the current directory are used.
#' @param vars Forcing names to visualize. Defaults to
#'   `c("pre", "pet", "tmin", "tmax", "tavg")`.
#' @examples
#' visualize_meteo_mosaic(domains = c("domain_1", "domain_2"))
#' @export
visualize_meteo_mosaic <- function(domains = NULL,
                                   vars = c("pre", "pet", "tmin", "tmax", "tavg")) {
  library(terra)
  library(jsonlite)

  if (is.null(domains)) {
    domains <- dir(pattern = "domain", full.names = TRUE)
  }
  if (length(domains) == 0) {
    stop("No domain folders found")
  }

  fun_map <- c(pre = "sum", pet = "sum", tmin = "mean", tmax = "mean", tavg = "mean")
  mosaics <- list()

  for (v in vars) {
    message("Processing variable: ", v)
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
      r_month <- daily_to_monthly(r, fun = fun_map[[v]])
      r_ann <- annual_mean(r_month, fun = fun_map[[v]])
      rasters[[length(rasters) + 1]] <- r_ann
    }
    mosaics[[v]] <- do.call(mosaic, rasters)
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 5))
  for (v in vars) {
    plot(mosaics[[v]], main = v)
  }
  if (length(vars) < 6) {
    for (i in seq_len(6 - length(vars))) plot.new()
  }

  invisible(mosaics)
}
