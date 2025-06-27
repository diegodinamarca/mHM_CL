library(terra)
library(tidyverse)
library(stars)
library(sf)
library(jsonlite)
source("R/utils.r")

# Read and merge tables with simulated and observed runoff
df = tibble()
for (i in 1:7){
  domain_path = paste0("../domain_zone_",i)
  x = read_csv(file.path(domain_path, "out", "streamflow","streamflow_data.csv"))
  df = bind_rows(df, x)
}

# gauge stations with unrealistic simulated streamflow
fault.ID = df %>% filter(Q_sim_m3s > 10^6 | Q_sim_m3s < 0) %>% 
  pull(ID) %>% 
  unique
fault.ID
# filter problematic stations
df = df %>% filter(!(ID %in% fault.ID))

# set inputs for metrics calculation
start_date = as.Date("1980-01-01")
end_date = as.Date("2020-12-31")
out_folder = "../domain_Chile/out/metrics"
dir.create(out_folder)
av_threshold = 80
config_path = "preprocess_config.json"
out.figs = "../FIGS"

# which map use as background
map.file = "../domain_Chile/out/Q.nc"

# uses simulated streamflow (with routing) to calculate metrics for each gauge
metrics <- evaluate_station_metrics(
  df = df,
  config_path,
  start_date,
  end_date,
  out_folder,
  av_threshold
)
# metrics = metrics %>% filter(kge > -20)

# Convert to sf
metrics_sp <- metrics %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326)


# calculate annual mean
r = rast(map.file) %>% annual_mean(fun = "sum")

# plot R2
p1 = plot_metric_map(r, metrics_sp, metric = "r2", 
                     # raster_palette = "magma",
                     # point_palette = "cividis",
                     raster_name = "Q [mm]",
                     metric_name = bquote(R^2)
)+
  # ggtitle(bquote("Mean annual accummulated total runoff and "*R^2*" of daily streamflow at gauges"))+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5))

# plot KGE
p2 = plot_metric_map(r, metrics_sp, metric = "kge", 
                     # raster_palette = "magma",
                     # point_palette = "cividis",
                     raster_name = "Q [mm]",
                     metric_name = "KGE"
)+
  # ggtitle(bquote("Mean annual accummulated total runoff and "*R^2*" of daily streamflow at gauges"))+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5))

p.mat = ggpubr::ggarrange(p1,p2, nrow = 1, labels = c("a","b"))
p.mat
ggsave(filename = file.path(out.figs,"gauge_metrics.png"), width = 10, height = 15, plot = p.mat)

# Cummulative probability distribution plot
# Choose the metric to plot
metric_name <- "kge"

# Prepare the data
plot_data <- metrics %>%
  ungroup() %>% 
  select(ID, !!sym(metric_name)) %>%
  arrange(desc(!!sym(metric_name))) %>%
  mutate(
    rank = row_number(),
    pct_above = 100 * (rank - 1) / n()
  )

# Plot
p.cdf = ggplot(plot_data, aes(x = !!sym(metric_name), y = 100 - pct_above)) +
  geom_line(color = "blue") +
  # geom_point(size = 1, alpha = 0.7) +
  labs(
    x = metric_name,
    y = "% of IDs with metric ≥ x",
    title = paste("Cumulative distribution of", metric_name)
  ) +
  theme_bw()
p.cdf
ggsave(file.path(out.figs, "cummulative_distribution_KGE_gauges.png"), width = 4, height = 4, plot = p.cdf)

# Correlation of streamflow with and without rounting

# labels for facet
labels = c("m3s" = "Runoff with routing (m3/s)",
           "mm" = "Runoff w/o routing (mm)")

# Daily streamflow plot
p1 = df %>% filter(!is.na(Q_obs_m3s)) %>% 
  pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var","source","unit")) %>% 
  mutate(var = str_c(var, "_", source)) %>% 
  select(ID, date, var, unit, value) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  ggplot(aes(x = Q_obs, y = Q_sim))+
  facet_wrap(.~unit, ncol = 2, scales = "free", 
             labeller = as_labeller(labels))+
  geom_point()+
  ggpubr::stat_cor()+
  # geom_abline(color = "grey")+
  ggtitle("Daily streamflow")

# calculate monthly streamflow
df.month = df %>% filter(!is.na(Q_obs_m3s)) %>% 
  mutate(date = floor_date(date, unit = "months")) %>% 
  group_by(ID, date) %>%
  summarise(count = n(),
            Q_obs_mm = sum(Q_obs_mm),
            Q_sim_mm = sum(Q_sim_mm),
            Q_obs_m3s = mean(Q_obs_m3s),
            Q_sim_m3s = mean(Q_sim_m3s))

# monthly streamflow plot
p2 = df.month %>% 
  filter(count > 28) %>% # Elimina las estaciones con menos de 25 observaciones al mes
  pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var","source","unit")) %>% 
  mutate(var = str_c(var, "_", source)) %>% 
  select(ID, date, var, unit, value) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  ggplot(aes(x = Q_obs, y = Q_sim))+
  facet_wrap(.~unit, ncol = 2, scales = "free",
             labeller = as_labeller(labels))+
  geom_point()+
  ggpubr::stat_cor()+
  # geom_smooth(color = "blue")+
  # geom_abline(color = "grey")+
  ggtitle("Monthly streamflow")

# plot grid
p.streamflow = ggpubr::ggarrange(p1,p2, ncol = 1, labels = c("a","b"))
p.streamflow
ggsave(filename = file.path(out.figs,"streamflow with and without routing.png"), height = 8, width = 8)


# cuencas pluviales y nivales
domain_path = "../domain_zone_7"
camels = read_sf("../DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/camels_cl_boundaries/camels_cl_boundaries.shp")
basins = c(9414001,7339001,5748001)
basins = c(9414001,7339001,5710001)

roi = camels %>% filter(gauge_id %in% basins)

# === Load config and paths ===
config_path <- file.path(domain_path, "preprocess_config.json")
config <- read_json(config_path)

streamflow_data_file <- config$streamflow_data_file_mm
gauges_file <- config$fluv_station_file

# === Read spatial data ===
roi = camels %>% filter(gauge_id %in% basins)
gauges <- read_csv(gauges_file) %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326)

# === Filter gauges by ROI if needed ===
sf::sf_use_s2(FALSE)
inter = st_intersection(roi, gauges) %>% as_tibble %>% select(gauge_id, gauge_name, ID)

df.basins = df %>% left_join(inter, by = "ID") %>% 
  filter(!is.na(gauge_id)) %>% 
  filter(!is.na(Q_obs_m3s))

df.basins = df.basins %>% filter(gauge_id == ID)

p1 = df.basins %>% 
  pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var","source","unit")) %>% 
  mutate(var = str_c(var, "_", unit)) %>% 
  select(gauge_id, gauge_name, ID, date, var, source, value) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  filter(source == "sim") %>% 
  ggplot(aes(x = Q_mm, y = Q_m3s))+
  facet_wrap(.~gauge_name, ncol = 3, scales = "free")+
  geom_point()+
  ggpubr::stat_cor()+
  labs(y = "Q with routing (m3/s)", x = "Q without routing (mm)")+
  # geom_abline(color = "grey")+
  ggtitle("Daily streamflow")


# calculate monthly streamflow
df.basins.month = df.basins %>% filter(!is.na(Q_obs_m3s)) %>% 
  mutate(date = floor_date(date, unit = "months")) %>% 
  group_by(gauge_id, gauge_name, ID, date) %>%
  summarise(count = n(),
            Q_obs_mm = sum(Q_obs_mm),
            Q_sim_mm = sum(Q_sim_mm),
            Q_obs_m3s = mean(Q_obs_m3s),
            Q_sim_m3s = mean(Q_sim_m3s))

# monthly streamflow plot
p2 = df.basins.month %>% 
  filter(count > 28) %>% # Elimina las estaciones con menos de 25 observaciones al mes
  pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var","source","unit")) %>% 
  mutate(var = str_c(var, "_", unit)) %>% 
  select(gauge_id, gauge_name, ID, date, var, source, value) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  filter(source == "sim") %>% 
  ggplot(aes(x = Q_mm, y = Q_m3s))+
  facet_wrap(.~gauge_name, ncol = 3, scales = "free")+
  geom_point()+
  ggpubr::stat_cor()+
  labs(y = "Q with routing (m3/s)", x = "Q without routing (mm)")+
  # geom_abline(color = "grey")+
  ggtitle("Monthly streamflow")


p.streamflow = ggpubr::ggarrange(p1,p2, ncol = 1, labels = c("a","b"))
p.streamflow
ggsave(filename = file.path(out.figs,"simulated streamflow routing by basin.png"), height = 8, width = 12)

