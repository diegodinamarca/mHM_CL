library(terra)
library(tidyverse)
library(stars)
library(sf)
library(yaml)
library(ggpubr)
library(lubridate)
source("R/utils.r")

# ============================
# FUNCIONES
# ============================

read_and_merge_streamflow <- function(config_name = "preprocess_config.yaml") {
  df <- tibble()
  # browser()
  if (use_multiple_domains) {
    dir.list <- dir("..", pattern = "domain_zone", full.names = TRUE)
    for (domain in dir.list) {
      config_path <- file.path(domain, config_name)
      config <- read_yaml(config_path)
      x <- read_csv(file.path(domain, config$out_folder, "streamflow", "streamflow_data.csv"))
      df <- bind_rows(df, x)
    }
  } else {
    config_path <- file.path(domain_path, config_name)
    config <- read_yaml(config_path)
    df <- read_csv(file.path(domain_path, config$out_folder, "streamflow", "streamflow_data.csv"))
  }
  
  fault.ID <- df %>% filter(Q_sim_m3s > 1e6 | Q_sim_m3s < 0) %>% pull(ID) %>% unique()
  df <- df %>% filter(!(ID %in% fault.ID))
  
  return(list(df = df, config = config))
}

setup_output_paths <- function(config) {
  if (use_multiple_domains) {
    out_figs <- "../FIGS"
    out_folder <- "../domain_Chile/out/metrics"
    map_file <- "../domain_Chile/out/Q.nc"
  } else {
    out_figs <- file.path(domain_path, "FIGS")
    out_folder <- file.path(domain_path, config$out_folder, "metrics")
    map_file <- file.path(domain_path, config$out_folder, "Q", "Q_month.tif")
  }
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)
  return(list(out_figs = out_figs, out_folder = out_folder, map_file = map_file))
}

calculate_metrics <- function(df, config, out_folder) {
  metrics <- evaluate_station_metrics(
    df = df,
    config_path = "preprocess_config.yaml",
    start_date = start_date,
    end_date = end_date,
    out_folder = out_folder,
    av_threshold = av_threshold
  )
  return(metrics)
}

plot_metrics_maps <- function(metrics, map_file, out_figs, 
                              plot_name = "gauge_metrics.png",
                              r2_threshold = NULL, 
                              kge_threshold = NULL) {
  # Convert metrics to sf
  # browser()
  dates_raw = str_sub(names(rast(map_file)), start = 2) 
  dates_vec = as.Date(paste0(dates_raw, ".01"), format = "%Y.%m.%d")
  metrics_sp <- metrics %>% st_as_sf(coords = c("LON", "LAT"), crs = 4326)
  r <- rast(map_file) %>% annual_mean(fun = "sum", dates_vector = dates_vec)
  
  # Filtrar según thresholds si están definidos
  metrics_r2 <- metrics_sp
  metrics_kge <- metrics_sp
  if (!is.null(r2_threshold)) {
    metrics_r2 <- metrics_r2 %>% filter(r2 >= r2_threshold)
  }
  if (!is.null(kge_threshold)) {
    metrics_kge <- metrics_kge %>% filter(kge >= kge_threshold)
  }
  
  # Gráfico R2
  p1 <- plot_metric_map(r, metrics_r2, metric = "r2", raster_name = "Q [mm]", metric_name = bquote(R^2)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Gráfico KGE
  p2 <- plot_metric_map(r, metrics_kge, metric = "kge", raster_name = "Q [mm]", metric_name = "KGE") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Plot conjunto
  p.mat <- ggarrange(p1, p2, nrow = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.mat, width = 10, height = 5)
}

plot_metrics_boxplot <- function(metrics, metric_name, out_figs, plot_name = NULL) {
  if (is.null(plot_name)) {
    plot_name <- paste0("boxplot_", metric_name, ".png")
  }
  
  p <- ggplot(metrics, aes(y = .data[[metric_name]])) +
    geom_boxplot(fill = "lightblue", outlier.shape = 21, outlier.fill = "red") +
    labs(y = metric_name, x = NULL, title = paste("Boxplot of", metric_name))
  
  ggsave(file.path(out_figs, plot_name), plot = p, width = 4, height = 6)
}


plot_cdf <- function(metrics, metric_name, out_figs, plot_name = NULL) {
  plot_data <- metrics %>%
    ungroup() %>%
    select(ID, !!sym(metric_name)) %>%
    arrange(desc(!!sym(metric_name))) %>%
    mutate(rank = row_number(), pct_above = 100 * (rank - 1) / n())
  
  p.cdf <- ggplot(plot_data, aes(x = !!sym(metric_name), y = 100 - pct_above)) +
    geom_line(color = "blue") +
    labs(x = metric_name, y = "% of IDs with metric ≥ x", title = paste("Cumulative distribution of", metric_name)) +
    theme_bw()
  
  if (is.null(plot_name)) {
    plot_name <- paste0("cummulative_distribution_", metric_name, "_gauges.png")
  }
  
  ggsave(file.path(out_figs, plot_name), plot = p.cdf, width = 4, height = 4)
}

plot_streamflow_comparison <- function(df, out_figs, plot_name = "streamflow_with_and_without_routing.png") {
  labels <- c("m3s" = "Runoff with routing (m3/s)", "mm" = "Runoff w/o routing (mm)")
  
  p1 <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", source)) %>%
    select(ID, date, var, unit, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    ggplot(aes(x = Q_obs, y = Q_sim)) +
    facet_wrap(. ~ unit, ncol = 2, scales = "free", labeller = as_labeller(labels)) +
    geom_point() +
    ggpubr::stat_cor() +
    ggtitle("Daily streamflow")
  
  df.month <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(ID, date) %>%
    summarise(count = n(),
              Q_obs_mm = sum(Q_obs_mm),
              Q_sim_mm = sum(Q_sim_mm),
              Q_obs_m3s = mean(Q_obs_m3s),
              Q_sim_m3s = mean(Q_sim_m3s))
  
  p2 <- df.month %>%
    filter(count > 28) %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", source)) %>%
    select(ID, date, var, unit, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    ggplot(aes(x = Q_obs, y = Q_sim)) +
    facet_wrap(. ~ unit, ncol = 2, scales = "free", labeller = as_labeller(labels)) +
    geom_point() +
    ggpubr::stat_cor() +
    ggtitle("Monthly streamflow")
  
  p.streamflow <- ggpubr::ggarrange(p1, p2, ncol = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.streamflow, width = 8, height = 8)
}

plot_basin_streamflow_comparison <- function(df, config, out_figs, basins, single_gauge = TRUE, plot_name = "simulated_streamflow_routing_by_basin.png") {
  # browser()
  camels <- read_sf("../DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/camels_cl_boundaries/camels_cl_boundaries.shp")
  gauges_file <- config$fluv_station_file
  gauges <- read_sf(gauges_file)
  
  roi <- camels %>% filter(gauge_id %in% basins)
  sf::sf_use_s2(FALSE)
  inter <- st_intersection(roi, gauges) %>% 
    as_tibble() %>% 
    select(gauge_id, gauge_name, ID)
  
  df.basins <- df %>% 
    left_join(inter, by = "ID") %>% 
    filter(!is.na(gauge_id), !is.na(Q_obs_m3s))
  
  if (single_gauge) {
    df.basins <- df.basins %>% filter(gauge_id == ID)
  }
  area = read_domain_area(domain_path, config_name) 
  df.basins = df.basins %>%  mutate(Q_sim_m3s = 60*60*24*Q_sim_m3s/(area*1000))
  p1 <- df.basins %>% 
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>% 
    mutate(var = str_c(var, "_", unit)) %>% 
    select(gauge_id, gauge_name, ID, date, var, source, value) %>% 
    pivot_wider(names_from = var, values_from = value) %>% 
    filter(source == "sim") %>% 
    ggplot(aes(x = Q_mm, y = Q_m3s)) +
    facet_wrap(. ~ gauge_name, ncol = 3, scales = "free") +
    geom_point() +
    ggpubr::stat_cor() +
    labs(y = "Q with routing [mm]", x = "Q without routing [mm]") +
    ggtitle("Daily streamflow")
  
  df.basins.month <- df.basins %>% 
    mutate(date = floor_date(date, "month")) %>% 
    group_by(gauge_id, gauge_name, ID, date) %>%
    summarise(count = n(),
              Q_obs_mm = sum(Q_obs_mm),
              Q_sim_mm = sum(Q_sim_mm),
              Q_obs_m3s = mean(Q_obs_m3s),
              Q_sim_m3s = sum(Q_sim_m3s),
              .groups = "drop")
  
  p2 <- df.basins.month %>% 
    filter(count > 28) %>% 
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>% 
    mutate(var = str_c(var, "_", unit)) %>% 
    select(gauge_id, gauge_name, ID, date, var, source, value) %>% 
    pivot_wider(names_from = var, values_from = value) %>% 
    filter(source == "sim") %>% 
    ggplot(aes(x = Q_mm, y = Q_m3s)) +
    facet_wrap(. ~ gauge_name, ncol = 3, scales = "free") +
    geom_point() +
    ggpubr::stat_cor() +
    labs(y = "Q with routing [mm]", x = "Q without routing [mm]") +
    ggtitle("Monthly streamflow")
  
  p.streamflow <- ggpubr::ggarrange(p1, p2, ncol = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.streamflow, width = 12, height = 8)
}

# ============================
# CONFIGURACIÓN GENERAL
# ============================

use_multiple_domains <- FALSE
domain_path <- "../domain_7339001"
start_date <- as.Date("1980-01-01")
end_date <- as.Date("2020-12-31")
av_threshold <- 50
suffix = "_default"
config_name = "preprocess_config_default.yaml"

# ============================
# WORKFLOW PRINCIPAL
# ============================

main <- function() {
  library(ggpubr)
  # browser()
  res <- read_and_merge_streamflow(config_name)
  df <- res$df
  config <- res$config
  # area = read_domain_area(domain_path, config_name)
  # df %>% mutate(qsim_route = 60*60*24*Q_sim_m3s/(area*1000))
  paths <- setup_output_paths(config)
  out_figs <- paths$out_figs
  out_folder <- paths$out_folder
  map_file <- paths$map_file
  
  metrics <- calculate_metrics(df, config, out_folder)
  metrics
  # plot_metrics_maps(metrics, map_file, out_figs, 
  #                   plot_name = paste0("only_valid_gauge_metrics",suffix,".png"),
  #                   r2_threshold = 0,
  #                   kge_threshold = 0)
  plot_metrics_maps(metrics, map_file, out_figs, 
                    plot_name = paste0("gauge_metrics",suffix,".png"))
  plot_metrics_boxplot(metrics, "r2", out_figs, plot_name = paste0("boxplot_r2",suffix,".png"))
  plot_metrics_boxplot(metrics, "kge", out_figs, plot_name = paste0("boxplot_kge",suffix,".png"))
  
  plot_cdf(metrics, "kge", out_figs, plot_name = paste0("cdf_kge",suffix,".png"))
  plot_streamflow_comparison(df, out_figs, plot_name = paste0("streamflow_comparison",suffix,".png"))
  
  basins <- c(9414001, 11337001, 7339001, 5710001)
  plot_basin_streamflow_comparison(df, config, out_figs, basins, single_gauge = TRUE, 
                                   plot_name = paste0("streamflow_basin",suffix,".png"))
}

main()

