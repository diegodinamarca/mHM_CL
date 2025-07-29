library(terra)
library(tidyverse)
library(stars)
library(sf)
library(yaml)
library(ggpubr)
library(lubridate)
source("R/utils.r")

# ---------------------------------------------------------------------------
# User configurable parameters
# ---------------------------------------------------------------------------
# domain_path       : Directory containing the mHM domain outputs
# config_name       : YAML configuration file produced during preprocessing
# use_multiple_domains : Set to TRUE to merge results from several domains
# start_date/end_date  : Date range to compute streamflow performance metrics
# av_threshold      : Minimum percentage of available observations required
#                     to include a gauge in the evaluation
# suffix            : String appended to the names of generated figures
# out.img           : Path for the summary image of annual means
# ---------------------------------------------------------------------------
domain_path = "../domain_zone_5"
config_name = "preprocess_config_default.yaml"
use_multiple_domains <- FALSE
start_date <- as.Date("1980-01-01")
end_date <- as.Date("2020-12-31")
av_threshold <- 50
suffix <- "_default"
out.img = file.path(domain_path, "FIGS/annual_output_default.png")
dir.create(file.path(domain_path, "FIGS"))

# extraer y escribir un archivo nc de una variable
variables = c("snowpack","SM_Lall","satSTW","aET","Q",
              "SWC_L01","SWC_L02","SWC_L03","SWC_L04","SWC_L05","SWC_L06")
for (j in c(6:11)){
  print(variables[j])
  write_output(domain_path, var_name = variables[j], ts = "month", roi_mask = TRUE, config_name = config_name)
}

var.clim = c("pre","pet")
for (j in c(1:2)){
  print(var.clim[j])
  write_clim(domain_path, var_name = var.clim[j], ts = "month", roi_mask = TRUE, config_name = config_name)
}

# visualizar todos los outputs y forcings de un dominio
calculate_annual_mean(domain_path, mask_roi = TRUE, config_name = config_name, process.clim = TRUE)
visualize_annual_mean(domain_path, mask_roi = TRUE, config_name = config_name, png_filename = out.img)

# Extract gridded runoff values (mm) and routing generated runoff (m3/s)
config_path <- file.path(domain_path, config_name)
config = read_yaml(config_path)
df.mm <- get_qmm_table(domain_path, crop_to_roi = TRUE, config_name = config_name)

discharge = file.path(domain_path, config$out_folder, "discharge.nc")
df <- get_qm3s_table(discharge)

gauges_time_range = df %>% group_by(ID) %>% 
  filter(!is.na(Q_obs)) %>% 
  filter(date >= "1980-01-01") %>%
  summarise(start = min(date),
            end = max(date)) %>% 
  mutate(years = (as.numeric(end-start)/365)) %>% 
  # ungroup() %>% 
  arrange(desc(years))
print(gauges_time_range)

df.full = df.mm %>% full_join(df, by = c("ID","date"), suffix = c("_mm","_m3s"))
dir.create(file.path(domain_path, config$out_folder,"streamflow"))
df.full %>% write_csv(file.path(domain_path, config$out_folder, "streamflow","streamflow_data.csv"))


# --------------------------------------------------
# Visualise streamflow performance metrics
# --------------------------------------------------
res <- read_and_merge_streamflow(config_name)
df.stream <- res$df
config <- res$config
paths <- setup_output_paths(config)

metrics <- calculate_metrics(df.stream, config, paths$out_folder)
plot_metrics_maps(metrics, paths$map_file, paths$out_figs,
                  plot_name = paste0("gauge_metrics", suffix, ".png"))
plot_metrics_boxplot(metrics, "r2", paths$out_figs,
                     plot_name = paste0("boxplot_r2", suffix, ".png"))
plot_metrics_boxplot(metrics, "kge", paths$out_figs,
                     plot_name = paste0("boxplot_kge", suffix, ".png"))
plot_cdf(metrics, "kge", paths$out_figs,
         plot_name = paste0("cdf_kge", suffix, ".png"))
plot_streamflow_comparison(df.stream, paths$out_figs,
                           plot_name = paste0("streamflow_comparison", suffix, ".png"))

basins <- c(9414001, 11337001, 7339001, 5710001)
plot_basin_streamflow_comparison(df.stream, config, paths$out_figs, basins, single_gauge = TRUE,
                                 plot_name = paste0("streamflow_basin", suffix, ".png"))

