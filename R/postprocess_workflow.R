library(terra)
library(tidyverse)
library(stars)
library(sf)
library(yaml)
source("R/utils.r")

# set the path for the domain
domain_path = "../domain_zone_5"
config_name = "preprocess_config_default.yaml"
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

# # extract single variable from mhm output
# df.q = extract_roi_timeseries("domain_zone_4/OUT/mHM_Fluxes_States.nc", var_name = "Q", roi_file = basin)
# # extract all outputs from mhm output
# df = extract_roi_timeseries_all("domain_zone_4/OUT/mHM_Fluxes_States.nc", roi_file = basin)
# # process forcings to monthly, yearly and annual mean for visualization
# visualize_mosaic_full_outputs("domain_chile/OUT") # this function takes too long
