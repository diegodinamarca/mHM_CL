#' Coordinate selector UI
#'
#' Provides inputs to upload a CSV file with columns `lon` and `lat` or
#' enter a pair of coordinates manually. Selected points are displayed on a
#' `leaflet` map.
#'
#' @param id Module id.
#' @return A shiny UI definition.
#' @export
coordinate_selector_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(ns("file"), "Upload CSV", accept = ".csv"),
    shiny::numericInput(ns("lon"), "Longitude", value = NA),
    shiny::numericInput(ns("lat"), "Latitude", value = NA),
    leaflet::leafletOutput(ns("map"), height = "400px")
  )
}

#' Coordinate selector server
#'
#' Reads the uploaded CSV or the manually entered coordinates and returns
#' them as an `sf` object in EPSG:4326. Points are drawn on a `leaflet` map.
#'
#' @param id Module id.
#' @return A reactive expression with selected points (`sf` object).
#' @export
coordinate_selector_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # reactive sf object with selected points
    coords_sf <- reactive({
      if (!is.null(input$file)) {
        df <- tryCatch(
          utils::read.csv(input$file$datapath),
          error = function(e) NULL
        )
        if (!is.null(df) && all(c("lon", "lat") %in% names(df))) {
          sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
        } else {
          sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
        }
      } else if (!is.na(input$lon) && !is.na(input$lat)) {
        df <- data.frame(lon = input$lon, lat = input$lat)
        sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
      } else {
        sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
      }
    })

    output$map <- leaflet::renderLeaflet({
      leaflet::leaflet() %>%
        leaflet::addProviderTiles("OpenStreetMap")
    })

    observe({
      pts <- coords_sf()
      proxy <- leaflet::leafletProxy("map", session)
      proxy %>% leaflet::clearMarkers()
      if (nrow(pts) > 0) {
        xy <- sf::st_coordinates(pts)
        proxy %>% leaflet::addCircleMarkers(lng = xy[,1], lat = xy[,2])
        bb <- sf::st_bbox(pts)
        proxy %>% leaflet::fitBounds(bb$xmin, bb$ymin, bb$xmax, bb$ymax)
      }
    })

    return(coords_sf)
  })
}
