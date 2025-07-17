#' Timeseries plot UI
#'
#' @param id Module id
#' @return Shiny UI elements for the timeseries plot module
#' @export
#' @keywords internal
#'
#' @examples
#' timeseries_plot_ui("plot1")

timeseries_plot_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::dateRangeInput(
      ns("date_range"),
      label = "Select date range",
      start = NULL,
      end = NULL
    ),
    plotly::plotlyOutput(ns("timeseries_plot"))
  )
}

#' Timeseries plot server
#'
#' @param id Module id
#' @param timeseries_data A `data.frame` or reactive expression with columns `date` and `value`
#' @export
#' @return None
#'
#' @examples
#' timeseries_plot_server("plot1", reactive_df)

timeseries_plot_server <- function(id, timeseries_data) {
  moduleServer(id, function(input, output, session) {
    data_reactive <- if (shiny::is.reactive(timeseries_data)) timeseries_data else shiny::reactive(timeseries_data)

    # Update date range based on data
    shiny::observe({
      df <- data_reactive()
      shiny::req(df)
      if (!all(c("date", "value") %in% names(df))) {
        stop("Input data must contain 'date' and 'value' columns")
      }
      shiny::updateDateRangeInput(session, "date_range",
                                 start = min(df$date, na.rm = TRUE),
                                 end = max(df$date, na.rm = TRUE))
    }, priority = 1)

    filtered_data <- shiny::reactive({
      df <- data_reactive()
      shiny::req(df)
      start <- input$date_range[1]
      end <- input$date_range[2]
      if (!is.null(start) && !is.na(start) && !is.null(end) && !is.na(end)) {
        df <- dplyr::filter(df, date >= start & date <= end)
      }
      df
    })

    output$timeseries_plot <- plotly::renderPlotly({
      df <- filtered_data()
      shiny::req(nrow(df) > 0)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
        ggplot2::geom_line(color = "steelblue") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = NULL, y = NULL)
      plotly::ggplotly(p)
    })
  })
}

