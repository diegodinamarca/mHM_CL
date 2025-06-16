source("R/utils.r")
source("R/postprocess/postprocessing_utils.r")



for (i in 1:7){
  domain_path = paste0("../domain_zone_",i)
  process_meteo_variable(domain_path = domain_path, var_name = "pre", roi_mask = TRUE)
  process_meteo_variable(domain_path = domain_path, var_name = "pet", roi_mask = TRUE)
}

# set the path for the domain
domain_path = "../domain_zone_2"
# visualizar todos los outputs y forcings de un domain
visualize_annual_outputs(domain_path, mask_roi = TRUE)

# extraer y escribir un archivo nc de una variable
for (i in 1:7){
  domain_path = paste0("../domain_zone_",i)
  print(domain_path)
  for (j in 1:5){
    variables = c("snowpack","SM_Lall","satSTW","aET","Q")
    print(variables[j])
    write_output(domain_path, var_name = variables[j], ts = "month", roi_mask = TRUE)
  }
}
# lista de dominios (zonification of chile)
domain_folders = dir(pattern = "domain_zone", full.names = TRUE)
domain_folders
# Mosaic outputs across domains
mosaic_outputs(domain_folders)
# visualize mosaic outputs
# visualize_mosaic_outputs("domain_chile/OUT")

# custom roi
basin = "/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/DATA/SHP/cau_basin.geojson"
# visualize mosaic outputs on custom roi
visualize_mosaic_outputs("domain_chile/OUT", roi_file = basin)

# extract single variable from mhm output
df.q = extract_roi_timeseries("domain_zone_4/OUT/mHM_Fluxes_States.nc", +var_name = "Q", roi_file = basin)
# extract all outputs from mhm output
df = extract_roi_timeseries_all("domain_zone_4/OUT/mHM_Fluxes_States.nc", roi_file = basin)

# process forcings to monthly, yearly and annual mean for visualization
visualize_mosaic_full_outputs("domain_chile/OUT") # this function takes too long
