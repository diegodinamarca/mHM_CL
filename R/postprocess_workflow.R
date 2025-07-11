library(terra)
library(tidyverse)
library(stars)
library(sf)
library(yaml)
source("R/utils.r")


# process_pet_pr <- function(domain_path) {
#     process_meteo_variable(domain_path = domain_path, var_name = "pre", roi_mask = TRUE)
#     process_meteo_variable(domain_path = domain_path, var_name = "pet", roi_mask = TRUE)
# }
# process_pet_pr("9404001")

# for processing all the zones at one 
# for (i in 1:7){
#   domain_path = paste0("../domain_zone_",i)
#   process_pet_pr(domain_path)
# }

# for processing all basin ids at once
# domain_list = list_domain_ids()
# for (d in domein_list){
#   domain_path = d
#   process_pet_pr(domain_path)
# }

# set the path for the domain
domain_path = "../domain_zone_6"

# # visualizar todos los outputs y forcings de un dominio
# out.img = file.path(domain_path, "FIGS/annual_output_default.png")
# visualize_annual_outputs(domain_path, mask_roi = TRUE, filename = out.img)

# extraer y escribir un archivo nc de una variable
variables = c("snowpack","SM_Lall","satSTW","aET","Q",
              "SM_L01","SM_L02","SM_L03","SM_L04","SM_L05","SM_L06")
for (j in 5:5){
  print(variables[j])
  write_output(domain_path, var_name = variables[j], ts = "month", roi_mask = TRUE)
}

# # Procesar todas las zonas
# for (i in 1:7){
#   domain_path = paste0("../domain_zone_",i)
#   print(domain_path)
#   for (j in 1:11){
#     variables = c("snowpack","SM_Lall","satSTW","aET","Q",
#                   "SM_L01","SM_L02","SM_L03","SM_L04","SM_L05","SM_L06")
#     print(variables[j])
#     write_output(domain_path, var_name = variables[j], ts = "month", roi_mask = TRUE)
#   }
# }

# # lista de dominios
# domain_folders = dir(pattern = "domain_zone", full.names = TRUE)
# domain_folders
# # Mosaic outputs across domains
# mosaic_outputs(domain_folders)

# Extract gridded runoff values (mm) and routing generated runoff (m3/s)
config_path <- file.path(domain_path, "preprocess_config.yaml")
config = read_yaml(config_path)
df.mm <- get_qmm_table(domain_path, crop_to_roi = TRUE)

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

# # running for all the zones at once
# for (i in 7:7){
#   domain_path = paste0("../domain_zone_",i)
#   config_path <- file.path(domain_path, "preprocess_config.yaml")
#   config = read_yaml(config_path)
#   df.mm <- get_qmm_table(domain_path, crop_to_roi = TRUE)
#   discharge = file.path(domain_path, config$out_folder, "discharge.nc")
#   df <- get_qm3s_table(discharge)
#   df.full = df.mm %>% full_join(df, by = c("ID","date"), suffix = c("_mm","_m3s"))
#   df.full %>% filter(!is.na(Q_obs_mm))
#   dir.create(file.path(domain_path, config$out_folder,"streamflow"))
#   df.full %>% write_csv(file.path(domain_path, config$out_folder, "streamflow","streamflow_data.csv"))
# }

# visualization!
# for (i in 7:7){
#   # i=1
#   domain_path = paste0("../domain_zone_",i)
#   config_path <- file.path(domain_path, "preprocess_config.yaml")
#   config = read_yaml(config_path)
#   discharge = file.path(domain_path, config$out_folder, "discharge.nc")
#   df <- get_qm3s_table(discharge)
#   df %>% group_by(ID) %>% 
#     filter(!is.na(Q_obs)) %>% 
#     filter(date >= "1980-01-01") %>% 
#     summarise(start = min(date),
#               end = max(date)) %>% 
#     mutate(years = (as.numeric(end-start)/365)) %>% 
#     # ungroup() %>% 
#     arrange(desc(years)) %>% 
#     filter(ID == 2112006)
#   df %>% filter(!is.na(Q_obs)) %>% 
#     filter(ID == 1021001)
# }

# extract single variable from mhm output
df.q = extract_roi_timeseries("domain_zone_4/OUT/mHM_Fluxes_States.nc", var_name = "Q", roi_file = basin)
# extract all outputs from mhm output
df = extract_roi_timeseries_all("domain_zone_4/OUT/mHM_Fluxes_States.nc", roi_file = basin)
# process forcings to monthly, yearly and annual mean for visualization
visualize_mosaic_full_outputs("domain_chile/OUT") # this function takes too long
