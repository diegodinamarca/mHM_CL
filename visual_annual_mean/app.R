#' Shiny application integrating raster navigation, coordinate selection and
#' time series visualization modules.
#'
#' Allows users to navigate raster files, select coordinates either from a CSV
#' file or manual entry, extract a time series from a NetCDF file at the chosen
#' point and visualise it interactively.
#'
#' The application assumes all spatial data use EPSG:4326.
#'
#' Run with `shiny::runApp()`.

# Load required packages ---------------------------------------------------
library(shiny)
library(leaflet)
library(sf)
library(terra)
library(plotly)
library(ggplot2)
library(dplyr)

# Source modules and helper function --------------------------------------
source("visual_annual_mean/raster_navigation_module.R")
source("visual_annual_mean/coordinate_selector_module.R")
source("visual_annual_mean/extract_timeseries.R")
source("visual_annual_mean/timeseries_plot_module.R")

# Vector of raster files to visualise. Update the path/pattern as needed.
raster_files <- list.files("rasters", pattern = "\.(tif|nc)$", full.names = TRUE)

# Path to NetCDF file containing the time series information.
nc_file <- "mhm_fluxes_states.nc"

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Raster visualisation and time series extraction"),
  fluidRow(
    column(6,
           raster_navigation_ui("raster_nav")
    ),
    column(6,
           coordinate_selector_ui("coords")
    )
  ),
  hr(),
  timeseries_plot_ui("ts_plot")
)

# Server ------------------------------------------------------------------
server <- function(input, output, session) {

  # Initialise modules -----------------------------------------------------
  raster_react <- raster_navigation_server("raster_nav", raster_files)
  coords_react <- coordinate_selector_server("coords")

  # Overlay selected points on raster map ---------------------------------
  observe({
    pts <- coords_react()
    proxy <- leafletProxy("raster_nav-map")
    proxy %>% clearMarkers()
    if (nrow(pts) > 0) {
      xy <- st_coordinates(pts)
      proxy %>% addCircleMarkers(lng = xy[,1], lat = xy[,2], color = "red")
    }
  })

  # Extract time series at first selected point ---------------------------
  ts_data <- reactive({
    pts <- coords_react()
    req(nrow(pts) > 0)
    extract_timeseries_at_point(pts[1, ], nc_file)
  })

  timeseries_plot_server("ts_plot", ts_data)
}

# Run the application -----------------------------------------------------
shinyApp(ui, server)
