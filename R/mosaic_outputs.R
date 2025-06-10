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
