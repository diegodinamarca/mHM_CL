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
  browser()
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
  par(mfrow = c(2, 4), mar = c(4, 4, 2, 5))

  plot(pre_ann, main = "Precipitation", zlim = flux_range)
  plot(pet_ann, main = "Potential ET", zlim = flux_range)
  plot(aet_ann, main = "Actual ET", zlim = flux_range)
  plot(runoff_ann, main = "Runoff", zlim = flux_range)

  plot(sm_ann, main = "Soil moisture", zlim = state_range)
  plot(snow_ann, main = "Snowpack", zlim = state_range)
  plot(sat_ann, main = "GW level", zlim = state_range)
  plot.new()
}
