

#' Load annual mean forcings across multiple domains
#'
#' Reads daily meteorological NetCDF files from each domain, masks them with the
#' ROI specified in `preprocess_config.json`, aggregates to monthly and then to
#' long term annual mean, and finally mosaics the domains.
#'
#' @param domains Character vector with paths to the domain folders. Defaults to
#'   `dir(pattern = "domain_zone", full.names = TRUE)`.
#' @param vars Character vector of forcing names to load.
#' @return A named list of [`terra::SpatRaster`] objects containing the annual
#'   mean of each variable.
#' @noRd
load_meteo_mosaic_annual_means <- function(domains = dir(pattern = "domain_zone", full.names = TRUE),
                                           vars) {
  library(terra)
  library(jsonlite)
  
  if (length(domains) == 0) {
    stop("No domain folders found")
  }
  if (missing(vars)) {
    stop("'vars' must be provided")
  }
  
  fun_map <- c(pre = "sum", pet = "sum", tmin = "mean", tmax = "mean", tavg = "mean")
  out <- vector("list", length(vars))
  names(out) <- vars
  
  for (v in vars) {
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
    out[[v]] <- do.call(mosaic, rasters)
  }
  
  out
}

#' Load annual mean from mosaicked outputs
#'
#' Reads NetCDF files produced by `mosaic_outputs()` and returns the long term
#' annual mean for each variable. If a region of interest (ROI) is provided,
#' results are masked to that area.
#'
#' @param out_dir Directory containing the mosaicked NetCDF files.
#' @param vars Character vector of variables to load. When `NULL`, all `.nc`
#'   files in `out_dir` are used.
#' @param roi_file Optional path to a ROI file used to mask the outputs.
#' @return A named list of [`terra::SpatRaster`] objects with the annual mean of
#'   each variable.
#' @export
load_mosaic_annual_means <- function(out_dir, vars = NULL, roi_file = NULL) {
  library(terra)

  if (!dir.exists(out_dir)) {
    stop("Directory not found: ", out_dir)
  }

  if (is.null(vars)) {
    files <- list.files(out_dir, pattern = "\\.nc$", full.names = FALSE)
    vars <- sub("\\.nc$", "", files)
  }

  roi <- NULL
  if (!is.null(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  }

  out <- vector("list", length(vars))
  names(out) <- vars

  for (v in vars) {
    nc_file <- file.path(out_dir, paste0(v, ".nc"))
    if (!file.exists(nc_file)) {
      stop("NetCDF not found: ", nc_file)
    }
    r <- rast(nc_file)
    fun <- if (v %in% c("aET", "Q", "pre", "pet")) "sum" else "mean"
    ann <- annual_mean(r, fun)
    if (!is.null(roi)) {
      ann <- crop(ann, roi)
      ann <- mask(ann, roi)
    }
    out[[v]] <- ann
  }

  out
}

#' Visualize annual mean of mosaicked outputs
#'
#' Displays the annual mean of the selected variables in a multi panel plot.
#'
#' @inheritParams load_mosaic_annual_means
#' @return Invisibly returns the list of rasters used for plotting.
#' @examples
#' visualize_mosaic_outputs("domain_chile/OUT", roi_file = "roi.shp")
#' @export
visualize_mosaic_outputs <- function(out_dir, vars = NULL, roi_file = NULL) {
  rasters <- load_mosaic_annual_means(out_dir, vars, roi_file)
  vars <- names(rasters)
  n <- length(rasters)
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 5))

  for (v in vars) {
    plot(rasters[[v]], main = v)
  }

  invisible(rasters)
}

#' Visualize mosaicked model outputs and forcings
#'
#' Combines annual mean forcings from multiple domains with the mosaicked model
#' outputs produced by `mosaic_outputs()` and displays them using the same 2x5
#' layout as `visualize_annual_outputs()`.
#'
#' @param out_dir Directory containing the mosaicked model NetCDF files.
#' @param domains Character vector with paths to the domain folders used to load
#'   the meteorological forcings. Defaults to
#'   `dir(pattern = "domain_zone", full.names = TRUE)`.
#' @param roi_file Optional ROI used to mask the model outputs.
#' @param vars Character vector of model output variables to visualize. Defaults
#'   to `c("snowpack", "SM_Lall", "satSTW", "aET", "Q")`.
#' @return Invisibly returns the list of rasters used for plotting.
#' @examples
#' domain_folders <- dir(pattern = "domain_zone", full.names = TRUE)
#' visualize_mosaic_full_outputs(out_dir = "domain_chile/OUT",
#'                               domains = domain_folders)
#' @export
visualize_mosaic_full_outputs <- function(out_dir,
                                          domains = dir(pattern = "domain_zone", full.names = TRUE),
                                          roi_file = NULL,
                                          vars = c("snowpack", "SM_Lall", "satSTW", "aET", "Q")) {
  library(terra)

  mod <- load_mosaic_annual_means(out_dir, vars, roi_file)
  met <- load_meteo_mosaic_annual_means(domains, c("pre", "pet"))

  pre_ann <- met[["pre"]]
  pet_ann <- met[["pet"]]
  aet_ann <- mod[["aET"]]
  q_ann   <- mod[["Q"]]
  sm_ann  <- mod[["SM_Lall"]]
  snow_ann <- mod[["snowpack"]]
  sat_ann <- mod[["satSTW"]]

  disp_ann <- pre_ann - aet_ann

  flux_range <- range(c(values(pre_ann), values(pet_ann), values(aet_ann), values(q_ann)), na.rm = TRUE)
  state_range <- range(c(values(snow_ann), values(sm_ann), values(sat_ann)), na.rm = TRUE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(2, 5), mar = c(4, 4, 2, 5))

  plot(pre_ann, main = "Pre", zlim = flux_range)
  plot(pet_ann, main = "Pet", zlim = flux_range)
  plot(aet_ann, main = "aET", zlim = flux_range)
  plot(q_ann, main = "Q", zlim = flux_range)
  plot(disp_ann, main = "Pr-ET", zlim = flux_range)

  plot(sm_ann, main = "SM_Lall", zlim = state_range)
  plot(snow_ann, main = "snowpack", zlim = state_range)
  plot(sat_ann, main = "satSTW", zlim = state_range)
  plot.new()

  invisible(c(met, mod))
}
