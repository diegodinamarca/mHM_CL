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
