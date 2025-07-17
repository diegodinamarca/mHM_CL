#' Shiny application integrating raster navigation, coordinate selection and
#' time series visualization modules.
#'
#' Allows users to dynamically choose raster files and a NetCDF file,
#' navigate through the rasters, select coordinates either from a CSV file or
#' manual entry, extract a time series at the chosen point and visualise it
#' interactively.
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
source("raster_navigation_module.R")
source("coordinate_selector_module.R")
source("extract_timeseries.R")
source("timeseries_plot_module.R")


# UI ----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Raster visualisation and time series extraction"),
  fluidRow(
    column(6,
           fileInput("raster_files", "Select raster files",
                     multiple = TRUE,
                     accept = c(".tif", ".nc")),
           fileInput("nc_file", "Select NetCDF file", accept = ".nc")
    )
  ),
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

  raster_files_react <- reactive({
    if (is.null(input$raster_files)) character(0) else input$raster_files$datapath
  })

  nc_file_react <- reactive({
    req(input$nc_file)
    input$nc_file$datapath
  })

  # Initialise modules -----------------------------------------------------
  raster_react <- raster_navigation_server("raster_nav", raster_files_react)
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
    req(nrow(pts) > 0, nc_file_react())
    extract_timeseries_at_point(pts[1, ], nc_file_react())
  })

  timeseries_plot_server("ts_plot", ts_data)
}

# Run the application -----------------------------------------------------
shinyApp(ui, server)
