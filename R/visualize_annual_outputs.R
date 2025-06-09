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

#' Visualize annual mean model outputs and forcings
#'
#' Reads the standard mHM NetCDF outputs together with precipitation and
#' potential evapotranspiration forcings and displays the long term annual means
#' in a 2x4 panel plot.
#'
#' @param domain_path Path to the model domain containing the `OUT` and `meteo`
#'   folders.
#' @examples
#' visualize_annual_outputs("/path/to/domain")
#' @export
visualize_annual_outputs <- function(domain_path) {
  library(terra)

  out_file <- file.path(domain_path, "OUT", "mHM_Fluxes_States.nc")
  meteo_folder <- file.path(domain_path, "meteo")

  snowpack <- rast(out_file, subds = "snowpack")
  soilmoist <- rast(out_file, subds = "SM_Lall")
  satstw <- rast(out_file, subds = "satSTW")
  aet <- rast(out_file, subds = "aET")
  runoff <- rast(out_file, subds = "Q")

  pre <- rast(file.path(meteo_folder, "pre.nc"), subds = "pre")
  pet <- rast(file.path(meteo_folder, "pet.nc"), subds = "pet")

  pre_ann <- annual_mean(pre, "sum")
  pet_ann <- annual_mean(pet, "sum")
  aet_ann <- annual_mean(aet, "sum")
  runoff_ann <- annual_mean(runoff, "sum")

  snow_ann <- annual_mean(snowpack, "mean")
  sm_ann <- annual_mean(soilmoist, "mean")
  sat_ann <- annual_mean(satstw, "mean")

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
