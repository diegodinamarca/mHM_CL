#' Extract time series at a single point from a NetCDF file
#'
#' Returns the values of the first raster variable inside `nc_file` at the
#' location specified by `point`. The NetCDF must contain a temporal dimension
#' named `time`. If the time units are given as `days since YYYY-MM-DD` they are
#' converted to proper `Date` objects. The result can optionally be restricted to
#' a range defined by `start_date` and `end_date`.
#'
#' @param point [`sf::st_point`] in EPSG:4326 specifying the target coordinate.
#' @param nc_file Path to `mhm_fluxes_states.nc` or a similar NetCDF file.
#' @param start_date Optional start date to limit the returned series.
#' @param end_date Optional end date to limit the returned series.
#'
#' @return A data frame with columns `date` and `value`.
#' @export
extract_timeseries_at_point <- function(point, nc_file,
                                        start_date = NULL, end_date = NULL) {
  library(terra)
  library(sf)

  if (!file.exists(nc_file)) {
    stop("NetCDF file not found: ", nc_file)
  }

  if (inherits(point, "sf")) {
    geom <- sf::st_geometry(point)
  } else if (inherits(point, "sfc")) {
    geom <- point
  } else {
    stop("`point` must be an sf or sfc POINT object")
  }

  crs_pt <- sf::st_crs(geom)
  if (is.na(crs_pt)) {
    stop("`point` lacks CRS information")
  }
  if (crs_pt$epsg != 4326) {
    geom <- sf::st_transform(geom, 4326)
  }

  p_vect <- terra::vect(sf::st_as_sf(geom))
  r <- terra::rast(nc_file)
  time_vals <- terra::time(r)

  if (is.numeric(time_vals) && !is.null(attr(time_vals, "units"))) {
    units_attr <- attr(time_vals, "units")
    if (grepl("since", units_attr)) {
      origin <- sub(".*since\\s+", "", units_attr)
      time_vals <- as.Date(time_vals, origin = origin)
    } else {
      time_vals <- as.Date(time_vals, origin = "1970-01-01")
    }
  } else {
    time_vals <- as.Date(time_vals)
  }

  vals <- as.vector(terra::extract(r, p_vect, ID = FALSE))
  df <- data.frame(date = time_vals, value = vals)

  if (!is.null(start_date)) {
    start_date <- as.Date(start_date)
    df <- df[df$date >= start_date, ]
  }
  if (!is.null(end_date)) {
    end_date <- as.Date(end_date)
    df <- df[df$date <= end_date, ]
  }

  return(df)
}

