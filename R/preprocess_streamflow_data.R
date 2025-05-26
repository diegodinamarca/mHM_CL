library(jsonlite)
library(lubridate)
library(tidyverse)

# === Load configuration ===
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

streamflow_data_file <- config$streamflow_data_file
gauge_list <- config$gauge_list
gauge_folder <- config$gauge_folder

dir.create(gauge_folder, showWarnings = FALSE, recursive = TRUE)

# === Read CSV assuming the first column is date ===
df <- read_csv(streamflow_data_file)
colnames(df)[1] <- "Date"

for (i in seq_along(gauge_list)) {
  gauge_id <- gauge_list[[i]]
  if (!(gauge_id %in% colnames(df))) {
    cat(sprintf("Gauge ID %s not found in data.\n", gauge_id))
    next
  }
  
  series_valid <- df %>% select(Date, !!gauge_id) %>% drop_na()
  series_valid

  if (nrow(series_valid) == 0) {
    cat(sprintf("No valid data for gauge %s. Skipping.\n", gauge_id))
    next
  }
  
  start_date <- min(series_valid$Date)
  end_date <- max(series_valid$Date)
  
  
  full_range <- data.frame(Date = seq(start_date, end_date, by = "1 day"))
  filled_series <- full_range %>%
    left_join(series_valid, by = "Date") %>%
    mutate(Value = ifelse(is.na(.data[[gauge_id]]), -9999, .data[[gauge_id]])) %>%
    mutate(Hour = 0, Minute = 0)
  
  # Header lines
  output_lines <- c(
    sprintf("%s Gauge %d (daily discharge)", gauge_id, i),
    "nodata   -9999",
    "n       1       measurements per day [1, 1440]",
    sprintf("start  %s", format(start_date, "%Y %m %d %H %M")),
    sprintf("end    %s", format(end_date, "%Y %m %d %H %M"))
  )
  
  # Data lines
  data_lines <- filled_series %>%
    mutate(Line = sprintf("%4d  %02d  %02d  %02d  %02d   %8.3f",
                          year(Date), month(Date), day(Date), Hour, Minute, Value)) %>%
    pull(Line)
  
  output_lines <- c(output_lines, data_lines)
  
  # Write output
  output_file <- file.path(gauge_folder, paste0(gauge_id, ".txt"))
  writeLines(output_lines, output_file)
  
  cat(sprintf("Wrote file: %s\n", output_file))
}
