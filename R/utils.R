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

#' Extract ROI mean time series
#'
#' Computes the spatial mean of a variable inside a region of interest (ROI)
#' for every time step of `mHM_Fluxes_States.nc`.
#'
#' @param nc_path Path to `mHM_Fluxes_States.nc`.
#' @param var_name Name of the variable to extract.
#' @param roi_file Either a path to a vector file or a [`terra::SpatVector`]
#'   representing the ROI. If polygons are provided, the mean value inside each
#'   polygon is extracted. For points, the value at each point is returned.
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

  if (inherits(roi_file, "SpatVector")) {
    roi <- roi_file
  } else if (is.character(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else {
    stop("`roi_file` must be a path or a SpatVector object")
  }

  r <- rast(nc_path, subds = var_name)
  geom <- unique(geomtype(roi))
  if (length(geom) != 1) {
    stop("ROI must contain a single geometry type")
  }

  if (geom == "points") {
    vals <- extract(r, roi, ID = FALSE)
  } else if (geom %in% c("polygons", "lines")) {
    vals <- extract(r, roi, fun = mean, na.rm = TRUE, ID = FALSE)
  } else {
    stop("Unsupported geometry type: ", geom)
  }

  vals <- t(as.matrix(vals))
  tvec <- terra::time(r)
  df <- data.frame(time = tvec, vals)
  colnames(df)[-1] <- paste0("roi", seq_len(ncol(vals)))
  
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
#' @param roi_file Either a path to a vector file or a [`terra::SpatVector`]
#'   representing the ROI. If polygons are provided, the mean value inside each
#'   polygon is extracted. For points, the value at each point is returned.
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

  if (inherits(roi_file, "SpatVector")) {
    roi <- roi_file
  } else if (is.character(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else {
    stop("`roi_file` must be a path or a SpatVector object")
  }

  geom <- unique(geomtype(roi))
  if (length(geom) != 1) {
    stop("ROI must contain a single geometry type")
  }

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  vars <- names(nc$var)

  r <- rast(nc_path)
  time_vec <- as.Date(terra::time(r))
  res_list <- list()

  for (v in vars) {
    r <- rast(nc_path, subds = v)
    if (geom == "points") {
      vals <- extract(r, roi, ID = FALSE)
    } else if (geom %in% c("polygons", "lines")) {
      vals <- extract(r, roi, fun = mean, na.rm = TRUE, ID = FALSE)
    } else {
      stop("Unsupported geometry type: ", geom)
    }
    vals <- t(as.matrix(vals))
    colnames(vals) <- paste0(v, "_roi", seq_len(ncol(vals)))
    res_list[[v]] <- vals
  }

  res_mat <- do.call(cbind, res_list)
  df <- data.frame(time = time_vec, res_mat)
  
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
#'   `domain_path` or the file provided via `roi_file`.
#' @param roi_file Optional path to a ROI file. When supplied, the outputs are
#'   cropped and masked to this ROI instead of the one defined in the
#'   configuration file.
#' @examples
#' visualize_annual_outputs("/path/to/domain")
#' @export
visualize_annual_outputs <- function(domain_path, mask_roi = FALSE,
                                     roi_file = NULL){
  library(terra)
  library(jsonlite)
  roi <- NULL
  if (!is.null(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else if (mask_roi) {
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
    if (!is.null(roi_file)) {
      pre_ann <- crop(pre_ann, roi)
      pet_ann <- crop(pet_ann, roi)
      aet_ann <- crop(aet_ann, roi)
      runoff_ann <- crop(runoff_ann, roi)

      sm_ann <- crop(sm_ann, roi)
      snow_ann <- crop(snow_ann, roi)
      sat_ann <- crop(sat_ann, roi)
    }

    pre_ann <- mask(pre_ann, roi)
    pet_ann <- mask(pet_ann, roi)
    aet_ann <- mask(aet_ann, roi)
    runoff_ann <- mask(runoff_ann, roi)

    sm_ann <- mask(sm_ann, roi)
    snow_ann <- mask(snow_ann, roi)
    sat_ann <- mask(sat_ann, roi)
  }

  disp_ann <- pre_ann - aet_ann
  
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


get_scale_lim = function(r, round.digit = 0){
  values(r) %>%
    quantile(., c(0.1, 0.9), na.rm = TRUE) %>% 
    round(., digits = round.digit)
}


plot_raster_facet <- function(r, depth, var, limits, facet_col = "attributes") {
  color.scale = scale_fill_distiller(type = "seq", palette = 16, direction = 1,
                                     na.value = "transparent", 
                                     name = str_c(var),
                                     limits = limits,
                                     oob = scales::squish)
  stars_r = stars::st_as_stars(r)
  ggplot() +
    geom_stars(data = stars_r) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    facet_wrap(.~get(facet_col))+
    theme(legend.title = element_text(size = 12))
  
}

plot_raster <- function(r, depth, var, limits, palette = 16) {
  color.scale = scale_fill_distiller(type = "seq", palette = palette, direction = 1,
                                     na.value = "transparent", 
                                     name = str_c(var," ",depth),
                                     limits = limits,
                                     oob = scales::squish)
  ggplot() +
    geom_stars(data = stars::st_as_stars(r)) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    theme(legend.title = element_text(size = 12))
}

theme_grid = theme(
  panel.background = element_rect(fill = "snow", ),
  panel.grid = element_line(color = "grey"),
  plot.background = element_rect(fill = "white"))


