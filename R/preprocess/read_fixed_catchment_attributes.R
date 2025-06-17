x = "/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/catchment_attributes.csv"
x = "/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/catchment_attributes_fix.csv"

read.csv(x) %>% as_tibble
df = read.csv(x, sep = ";", skip = 1) %>% as_tibble() %>% 
  select(ID = gauge_id, gauge_name, LAT = outlet_camels_lat, LON = outlet_camels_lon)

df %>% write_csv("/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/DATA/TABLE/est_fluv/catchment_gauge_coords.csv")
