#' Raster navigation UI
#'
#' Creates the user interface for navigating through a list of raster files.
#'
#' @param id Module ID.
#'
#' @return A list of UI elements for inclusion in a Shiny app.
#' @export
raster_navigation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    leaflet::leafletOutput(ns("map")),
    shiny::div(
      style = "margin-top: 5px;",
      shiny::actionButton(ns("prev"), "Previous"),
      shiny::span(shiny::textOutput(ns("counter"), inline = TRUE),
                  style = "margin: 0 10px;"),
      shiny::actionButton(ns("next"), "Next")
    )
  )
}

#' Raster navigation server
#'
#' Server logic for navigating and displaying raster files on an interactive
#' map. The current raster and index are returned as reactive values for use in
#' other modules.
#'
#' @param id Module ID.
#' @param raster_files Character vector with paths to raster files (`.tif` or
#'   `.nc`).
#'
#' @return A list with reactive elements `raster` and `index`.
#' @export
raster_navigation_server <- function(id, raster_files) {
  moduleServer(id, function(input, output, session) {
    n <- length(raster_files)
    idx <- reactiveVal(1)

    observeEvent(input$prev, {
      if (idx() > 1) idx(idx() - 1)
    })
    observeEvent(input$next, {
      if (idx() < n) idx(idx() + 1)
    })

    current_raster <- reactive({
      r <- terra::rast(raster_files[idx()])
      terra::project(r, "EPSG:4326")
    })

    output$map <- leaflet::renderLeaflet({
      r <- current_raster()
      rng <- range(terra::values(r), na.rm = TRUE)
      pal <- leaflet::colorNumeric(viridisLite::viridis(256), rng,
                                   na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(r, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = rng)
    })

    output$counter <- shiny::renderText({
      sprintf("Raster %d de %d", idx(), n)
    })

    list(raster = current_raster, index = idx)
  })
}
