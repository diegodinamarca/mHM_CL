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
#' in the respective `preprocess_config.yaml` and merges all domains into
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
                           vars = c("snowpack","SM_Lall","satSTW","aET","Q",
                                    "SM_L01","SM_L02","SM_L03","SM_L04","SM_L05","SM_L06")) {
  library(terra)
  library(yaml)
  
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
      cfg_path <- file.path(d, "preprocess_config.yaml")
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- read_yaml(cfg_path)
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
#' `preprocess_config.yaml` is applied prior to mosaicking.
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
  library(yaml)
  
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
      cfg_path <- file.path(d, "preprocess_config.yaml")
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- read_yaml(cfg_path)
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
#'   If `NULL`, the `roi_file` defined in `preprocess_config.yaml` is used.
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
  library(yaml)
  
  out.format <- match.arg(out.format)
  ts <- match.arg(tolower(ts), c("month", "year"))
  
  config_path <- file.path(domain_path, "preprocess_config.yaml")
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  config <- read_yaml(config_path)
  
  nc_path <- file.path(domain_path, config$out_folder, "mHM_Fluxes_States.nc")
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }
  
  r <- rast(nc_path, subds = var_name)
  
  if (roi_mask) {
    if (is.null(roi_file)) {
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }
  
  if (ts == "month") {
    if (var_name %in% c("aET","Q")){
      r <- daily_to_monthly(r, fun = "sum")
    }else{
      r <- daily_to_monthly(r, fun = "mean")
    }
  } else if (ts == "year") {
    r <- monthly_to_yearly(r, fun = "mean")
  }
  
  out_dir <- file.path(domain_path, config$out_folder, var_name)
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
#' in a 2x5 panel plot.
#'
#' @param domain_path Path to the model domain containing the `OUT` and `meteo`
#'   folders.
#' @param mask_roi Logical. If `TRUE`, the outputs are masked using the
#'   `roi_file` specified in `preprocess_config.yaml` located inside
#'   `domain_path` or the file provided via `roi_file`.
#' @param roi_file Optional path to a ROI file. When supplied, the outputs are
#'   cropped and masked to this ROI instead of the one defined in the
#'   configuration file.
#' @param out_folder Folder where to save the output PNG figure. If `NULL`,
#'   defaults to "/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/FIGS".
#' @export
visualize_annual_outputs <- function(domain_path, mask_roi = FALSE,
                                     roi_file = NULL, filename = NULL, 
                                     out_folder = NULL, res_folder = NULL) {
  library(terra)
  library(yaml)
  
  # Internal function for annual mean
  annual_mean <- function(r, fun = "mean") {
    if (is.character(time(r))) {
      time(r) <- as.POSIXct(time(r), origin = "1970-01-01")
    }
    years <- format(time(r), "%Y")
    r_ann <- tapp(r, years, fun = fun, na.rm = TRUE)
    terra::mean(r_ann, na.rm = TRUE)
  }
  
  roi <- NULL
  if (!is.null(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else if (mask_roi) {
    config_path <- file.path(domain_path, "preprocess_config.yaml")
    if (!file.exists(config_path)) {
      stop("Configuration file not found: ", config_path)
    }
    config <- yaml::read_yaml(config_path)
    roi_path <- config$roi_file
    if (!grepl("^(/|[A-Za-z]:)", roi_path)) {
      roi_path <- file.path(domain_path, roi_path)
    }
    if (!file.exists(roi_path)) {
      stop("ROI file not found: ", roi_path)
    }
    roi <- vect(roi_path)
  }
  
  # Read datasets
  if (is.null(res_folder)){
    res_folder = config$out_folder
    res_file <- file.path(domain_path, res_folder, "mHM_Fluxes_States.nc")
  }else{
    res_file <- file.path(domain_path, res_folder, "mHM_Fluxes_States.nc")
  }
  meteo_folder <- file.path(domain_path, config$meteo_folder)
  
  cat("reading results in", res_file)
  
  snowpack <- rast(res_file, subds = "snowpack")
  soilmoist <- rast(res_file, subds = "SM_Lall")
  satstw <- rast(res_file, subds = "satSTW")
  aet <- rast(res_file, subds = "aET")
  runoff <- rast(res_file, subds = "Q")
  
  pre <- rast(file.path(meteo_folder, "pre.nc"), subds = "pre")
  pet <- rast(file.path(meteo_folder, "pet.nc"), subds = "pet")
  
  # Process layers
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
  
  # Crop and mask if ROI is provided
  if (!is.null(roi)) {
    pre_ann <- crop(pre_ann, roi) |> mask(roi)
    pet_ann <- crop(pet_ann, roi) |> mask(roi)
    aet_ann <- crop(aet_ann, roi) |> mask(roi)
    runoff_ann <- crop(runoff_ann, roi) |> mask(roi)
    sm_ann <- crop(sm_ann, roi) |> mask(roi)
    snow_ann <- crop(snow_ann, roi) |> mask(roi)
    sat_ann <- crop(sat_ann, roi) |> mask(roi)
  }
  
  disp_ann <- pre_ann - aet_ann
  
  flux_range <- range(c(values(pre_ann), values(pet_ann),
                        values(aet_ann), values(runoff_ann),
                        values(disp_ann)), na.rm = TRUE)
  
  state_range <- range(c(values(snow_ann), values(sm_ann), values(sat_ann)),
                       na.rm = TRUE)
  
  # Save PNG
  png_filename <- filename
  png(png_filename, width = 2000, height = 1000, res = 150)
  
  par(mfrow = c(2, 5), mar = c(4, 4, 2, 5))
  
  plot(pre_ann, main = "Precipitation", zlim = flux_range)
  plot(pet_ann, main = "Potential ET", zlim = flux_range)
  plot(aet_ann, main = "Actual ET", zlim = flux_range)
  plot(runoff_ann, main = "Runoff", zlim = flux_range)
  plot(disp_ann, main = "Pr - ET", zlim = flux_range)
  
  plot(sm_ann, main = "Soil Moisture", zlim = state_range)
  plot(snow_ann, main = "Snowpack", zlim = state_range)
  plot(sat_ann, main = "GW Level", zlim = state_range)
  plot.new()
  text(0.5, 0.5, "Summary Map", cex = 1.2)
  
  dev.off()
  message("Figure saved to: ", png_filename)
  
  # Plot to screen if interactive
  if (interactive()) {
    par(mfrow = c(2, 5), mar = c(4, 4, 2, 5))
    
    plot(pre_ann, main = "Precipitation", zlim = flux_range)
    plot(pet_ann, main = "Potential ET", zlim = flux_range)
    plot(aet_ann, main = "Actual ET", zlim = flux_range)
    plot(runoff_ann, main = "Runoff", zlim = flux_range)
    plot(disp_ann, main = "Pr - ET", zlim = flux_range)
    
    plot(sm_ann, main = "Soil Moisture", zlim = state_range)
    plot(snow_ann, main = "Snowpack", zlim = state_range)
    plot(sat_ann, main = "GW Level", zlim = state_range)
    plot.new()
    text(0.5, 0.5, "Summary Map", cex = 1.2)
  }
  
  invisible(list(
    pre = pre_ann,
    pet = pet_ann,
    aet = aet_ann,
    runoff = runoff_ann,
    disp = disp_ann,
    sm = sm_ann,
    snow = snow_ann,
    sat = sat_ann
  ))
}



#' Get plotting limits for a raster
#'
#' Computes quantile-based limits of the raster values which are typically used
#' to define the limits of a colour scale.
#' Computes the 10\% and 90\% quantiles of the raster values which are
#' typically used to define the limits of a colour scale.
#'
#' @param r A [`terra::SpatRaster`] object.
#' @param round.digit Integer indicating the number of digits to round the
#'   limits.
#' @param lower.lim Numeric lower quantile used when computing limits.
#' @param upper.lim Numeric upper quantile used when computing limits.
#'
#' @return Numeric vector of length two containing the lower and upper limits.
#' @keywords internal
get_scale_lim = function(r, round.digit = 0, lower.lim = 0.1, upper.lim = 0.9){
  values(r) %>%
    quantile(., c(lower.lim, upper.lim), na.rm = TRUE) %>%
    round(., digits = round.digit)
}

#' Plot raster facets
#'
#' Displays multiple layers of a raster in a faceted layout with a common
#' colour scale.
#'
#' @param r A [`terra::SpatRaster`] object.
#' @param var Character. Title for the legend.
#' @param limits Numeric vector providing the limits of the colour scale.
#' @param facet_col Attribute column used to create the facets.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_raster_facet <- function(r, var, limits, facet_col = "attributes") {
  library(stars)
  library(cowplot)
  color.scale = scale_fill_distiller(type = "seq", palette = 16, direction = 1,
                                     na.value = "transparent",
                                     name = var,
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

#' Plot a single raster
#'
#' Creates a map of one raster layer using a sequential palette.
#'
#' @param r A [`terra::SpatRaster`] to plot.
#' @param var Character string used as legend title.
#' @param limits Numeric vector with limits for the colour scale.
#' @param palette Integer palette number passed to `scale_fill_distiller`.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_raster <- function(r, var, limits, palette = 16) {
  color.scale = scale_fill_distiller(type = "seq", palette = palette, direction = 1,
                                     na.value = "transparent",
                                     name = var,
                                     limits = limits,
                                     oob = scales::squish)
  ggplot() +
    geom_stars(data = stars::st_as_stars(r)) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    theme(legend.title = element_text(size = 12))
}


#' Get streamflow table in millimetres. Extracts runoff values from gauge stations.
#'
#' Reads simulated runoff from `mHM_Fluxes_States.nc` together with CAMELS
#' observations and returns a table.
#'
#' @param domain_path Path to the model domain directory.
#' @param crop_to_roi Logical. If `TRUE`, only gauges inside the ROI defined in
#'   `preprocess_config.yaml` are considered.
#'
#' @return A tibble with columns `ID`, `date`, `Q_sim` and `Q_obs`.
#' @export
get_qmm_table <- function(domain_path, crop_to_roi = TRUE) {
  library(terra)
  library(sf)
  library(yaml)
  library(tidyverse)
  
  # === Load config and paths ===
  config_path <- file.path(domain_path, "preprocess_config.yaml")
  config <- read_yaml(config_path)
  
  streamflow_data_file <- config$streamflow_data_file_mm
  gauges_file <- config$fluv_station_file
  roi_file <- config$roi_file
  gauge_folder <- file.path(domain_path, config$gauge_folder)
  nc_path <- file.path(domain_path, config$out_folder, "mHM_Fluxes_States.nc")
  
  # === Read spatial data ===
  roi <- read_sf(roi_file)
  sf_use_s2(FALSE)
  gauges <- read_sf(gauges_file)
  
  # === Filter gauges by ROI if needed ===
  if (crop_to_roi) {
    gauge_list <- as.character(st_intersection(roi, gauges)$ID)
    message(length(gauge_list), " gauge stations found in ROI: ", paste(gauge_list, collapse = ", "))
  } else {
    gauge_list <- as.character(gauges$ID)
    message(length(gauge_list), " gauge stations used: ", paste(gauge_list, collapse = ", "))
  }
  
  gauges_roi <- gauges %>% filter(ID %in% gauge_list)
  
  # === Read Q simulation from NetCDF ===
  nc <- rast(nc_path, subds = "Q")
  names(nc) <- as.character(time(nc))  # Assign date names to layers
  
  val <- terra::extract(nc, gauges_roi)
  val$ID <- gauges_roi$ID
  
  df_sim <- val %>%
    as_tibble() %>%
    pivot_longer(cols = -ID, names_to = "date", values_to = "Q_sim") %>%
    mutate(date = as.Date(date))
  
  # === Read CAMELS observation ===
  df_obs <- read_csv(streamflow_data_file) %>%
    select(date, all_of(gauge_list)) %>%
    pivot_longer(cols = -date, names_to = "ID", values_to = "Q_obs") %>%
    mutate(ID = as.numeric(ID), date = as.Date(date))
  
  # === Join simulated and observed ===
  df <- full_join(df_sim, df_obs, by = c("ID", "date"))
  sf_use_s2(TRUE)
  
  return(df)
}

#' Get streamflow table in m\eqn{^3}/s. Reads discharge simulated with the routing
#' option of m
#'
#' Extracts simulated and observed discharge from a NetCDF file (discharge.nc)
#' generated by mHM
#'
#' @param nc_path Path to the NetCDF file containing variables `Qsim_<id>` and
#'   `Qobs_<id>`.
#'
#' @return A tibble with columns `ID`, `date`, `Q_sim` and `Q_obs`.
#' @export
get_qm3s_table <- function(nc_path) {
  library(tidyverse)
  library(ncdf4)
  
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))
  
  # Extract and convert time
  time <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  origin <- str_extract(time_units, "\\d{4}-\\d{2}-\\d{2}")
  dates <- as.Date(time / 24, origin = origin)  # assuming time is in hours
  
  # Get variable names
  var_names <- names(nc$var)
  sim_vars <- var_names[str_detect(var_names, "^Qsim_")]
  obs_vars <- var_names[str_detect(var_names, "^Qobs_")]
  
  # Extract gauge IDs
  sim_ids <- str_remove(sim_vars, "^Qsim_")
  obs_ids <- str_remove(obs_vars, "^Qobs_")
  common_ids <- intersect(sim_ids, obs_ids)
  
  # Function to extract time series
  extract_ts <- function(var) {
    vals <- ncvar_get(nc, var)
    vals[vals == -9999] <- NA
    vals
  }
  
  # Combine all into tibble
  result <- purrr::map_dfr(common_ids, function(id) {
    qsim <- extract_ts(paste0("Qsim_", id))
    qobs <- extract_ts(paste0("Qobs_", id))
    tibble(
      ID = as.numeric(id),
      date = dates,
      Q_sim = as.numeric(qsim),
      Q_obs = as.numeric(qobs)
    )
  })
  
  return(result)
}

#' Calculate hydrological error metrics
#'
#' Computes a range of performance metrics comparing observed and simulated
#' streamflow.
#'
#' @param obs Numeric vector of observations.
#' @param sim Numeric vector of simulations.
#' @param norm Character string passed to `hydroGOF::nrmse` indicating the
#'   normalisation method.
#'
#' @return A tibble with one row containing the calculated metrics.
#' @keywords internal
calculate_error_metrics <- function(obs, sim, norm = "sd") {
  # browser()
  N <- length(obs)
  bias <- mean(sim - obs)
  pbias <- hydroGOF::pbias(sim, obs)
  r2 <- hydroGOF::R2(sim, obs)
  rmse <- hydroGOF::rmse(sim, obs)
  ubrmse <- hydroGOF::ubRMSE(sim, obs)
  nrmse <- hydroGOF::nrmse(sim, obs, norm = norm)
  nse <-hydroGOF::NSE(sim, obs)
  kge_components <- hydroGOF::KGE(sim, obs, method="2009", out.type = "full")
  kge <- kge_components$KGE.value
  alpha <- kge_components$KGE.elements[["Alpha"]]
  beta <- kge_components$KGE.elements[["Beta"]]
  r <- kge_components$KGE.elements[["r"]]
  mae <- hydroGOF::mae(sim, obs)
  hbe <- (sum(sim)-sum(obs))/sum(obs)
  return(tibble(
    kge = kge,
    r2 = r2,
    bias = bias,
    pbias = pbias,
    rmse = rmse,
    ubrmse = ubrmse,
    nrmse = nrmse,
    nse = nse,
    alpha = alpha,
    beta = beta,
    hbe = hbe,
    r = r,
    mae = mae,
    N = N
  ))
}



#' Evaluate streamflow metrics at gauge stations
#'
#' Calculates performance metrics for each gauge station in the supplied data
#' frame and optionally writes the results to disk.
#'
#' @param df Data frame containing at least the columns `ID`, `date`,
#'   `Q_obs_mm`, `Q_sim_mm`, `Q_obs_m3s` and `Q_sim_m3s`.
#' @param config_path Path to `preprocess_config.yaml` with station metadata.
#' @param start_date Start date of the evaluation period.
#' @param end_date End date of the evaluation period.
#' @param out_folder Optional output folder where the table of metrics is
#'   written.
#' @param av_threshold Minimum percentage of available observations required to
#'   include a gauge.
#'
#' @return A tibble with the computed metrics and station coordinates.
#' @export
evaluate_station_metrics <- function(df,
                                     config_path = "preprocess_config.yaml",
                                     start_date = as.Date("1980-01-01"),
                                     end_date = as.Date("2020-12-31"),
                                     out_folder = NULL,
                                     av_threshold = 80) {
  # Read config file
  config <- yaml::read_yaml(config_path)
  # browser()
  calculate_non_na_percentage <- function(df, start_date = as.Date("1980-01-01"), end_date = as.Date("2020-12-31")) {
    df %>%
      filter(date >= start_date & date <= end_date) %>%
      group_by(ID) %>%
      summarise(
        total = n(),
        non_na_Q_obs_mm = sum(!is.na(Q_obs_mm)),
        non_na_Q_sim_mm = sum(!is.na(Q_sim_mm)),
        non_na_Q_obs_m3s = sum(!is.na(Q_obs_m3s)),
        non_na_Q_sim_m3s = sum(!is.na(Q_sim_m3s)),
        .groups = "drop"
      ) %>%
      mutate(
        perc_Q_obs_mm = 100 * non_na_Q_obs_mm / total,
        perc_Q_sim_mm = 100 * non_na_Q_sim_mm / total,
        perc_Q_obs_m3s = 100 * non_na_Q_obs_m3s / total,
        perc_Q_sim_m3s = 100 * non_na_Q_sim_m3s / total
      ) %>%
      select(ID, starts_with("perc_"))
  }
  # Calculate non-NA percentages (assumes the function is defined elsewhere)
  perc <- calculate_non_na_percentage(df)
  
  selected_gauges <- perc %>%
    filter(perc_Q_obs_m3s > av_threshold) %>%
    pull(ID)
  if (length(selected_gauges) == 0){
    selected_gauges = perc$ID
    message("All guages have less streamflow data than %", av_threshold, ", using all data available")
  }
  
  # Filter and compute metrics (assumes calculate_error_metrics() is defined)
  metrics <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    filter(date >= start_date, date <= end_date) %>%
    filter(ID %in% selected_gauges) %>%
    group_by(ID) %>%
    do(calculate_error_metrics(.$Q_obs_m3s, .$Q_sim_m3s))
  print(metrics)
  # Read station coordinates from config
  message("Reading config file. Fluvio stations coordinates read at: ", config_path)
  fluv_station <- read_sf(config$fluv_station_file)
  
  # Join coordinates
  metrics <- metrics %>%
    left_join(fluv_station, by = "ID")
  
  # Optionally write output
  if (!is.null(out_folder)) {
    dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
    write_csv(metrics, file.path(out_folder, "eval_metrics_stations.csv"))
  }
  
  return(metrics)
}

#' Plot raster background with station metrics
#'
#' Generates a map showing a raster layer together with station performance
#' metrics coloured by the chosen variable.
#'
#' @param rast A [`terra::SpatRaster`] used as background.
#' @param metrics_sf An [`sf`] object containing the metric values.
#' @param metric Name of the metric column to colour the points.
#' @param metric_name Optional title for the point legend.
#' @param raster_name Optional title for the raster legend.
#' @param raster_palette Viridis palette used for the raster fill.
#' @param point_palette Viridis palette used for the point colours.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_metric_map <- function(rast, metrics_sf, metric = "kge", metric_name = NULL,
                            raster_name = NULL,
                            raster_palette = "viridis", point_palette = "plasma") {
  # scale_fill_viridis_c() and scale_color_viridis_c()
  # support various palettes: "viridis", "magma", "plasma", "inferno", "cividis".
  # Convert raster to stars object
  stars_rast <- st_as_stars(rast)
  if (is.null(metric_name)){
    metric_name = metric
  }
  if (is.null(raster_name)){
    raster_name = varnames(rast)
  }
  # Base plot
  p <- ggplot() +
    geom_stars(data = stars_rast, na.rm = TRUE) +
    scale_fill_viridis_c(option = raster_palette, na.value = NA, name = raster_name) +
    geom_sf(data = metrics_sf, aes_string(color = metric), size = 2, na.rm = TRUE) +
    scale_color_viridis_c(option = point_palette, name = metric_name, na.value = NA) +
    # theme_minimal() +
    theme(legend.position = "right") +
    labs(x = "Longitude", y = "Latitude")
  
  return(p)
}

streamflow_data_range <- function(df, id = NULL) {
  if (!is.null(id))
    df %>% dplyr::select(date, all_of(id)) %>% 
    drop_na() %>% 
    summarise(n = n(),
              start = min(date),
              end = max(date),
              years = n/365)
}