# procesar DEM para mHM
{
  # HEADER --------------------------------------------
  #
  # Author: Diego Dinamarca
  # Email:  ddinamarcamuller@gmail.com
  # 
  # Date:
  #
  # Script Name: 
  #
  # Script Description: preparar datos de soilgrids para modelacion
  #
  # Notes:
  #
  #
  # INSTALL PACKAGES & LOAD LIBRARIES -----------------
  cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
  packages <- c("tidyverse",
                "terra", 
                "sf",
                "here",
                "jsonlite") # list of packages to load
  n_packages <- length(packages) # count how many packages are required
  
  new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
  
  # install missing packages
  if(length(new.pkg)){
    install.packages(new.pkg)
  }
  
  # load all requried libraries
  for(n in 1:n_packages){
    cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
    lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
    eval(parse(text = lib_load)) # evaluate the string to load the library
  }
  # SET WORKING DIRECTORY -----------------------------
  cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
  wd <- here::here()
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
  
  # SET OPTIONS ---------------------------------------
  cat("SETTING OPTIONS... \n\n", sep = "")
  options(scipen = 999) # turns off scientific notation
  options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
  
  
  # CONFLICTS ---------------------------------------------------------------
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("extract", "terra")
  conflicted::conflict_prefer("last", "dplyr")
  # LOAD FUNCTIONS ------------------------------------
}

# Prepare mask for chile
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")
morph_folder = config$morph_folder
lc_folder = config$lc_folder
lai_folder = config$lai_folder
roi_file = config$roi_file
basin.mask = FALSE


dem_file = file.path(morph_folder, "dem.asc")
geo_file = file.path(morph_folder, "geology_class.asc")
lc_file = file.path(lc_folder, "landcover.asc")
lai_file = file.path(lai_folder, "lai.nc")



dem = rast(dem_file)
geo = rast(geo_file)
lc = rast(lc_file)
lai = rast(lai_file)

dem.m = !is.na(dem)
geo.m = !is.na(geo)
lc.m = !is.na(lc)
lai.m = !is.na(lai[[1]])

roi_mask = dem.m + geo.m + lc.m + lai.m
max.val = max(values(roi_mask), na.rm = TRUE)
roi_mask[roi_mask != max.val] = NA
plot(roi_mask)
writeRaster(roi_mask, file.path(morph_folder, "roi_mask.tif"), overwrite = TRUE)


# apply mask to dem
roi_mask = rast(file.path(morph_folder, "roi_mask.tif"))
if (basin.mask){
  roi_mask = read_sf(roi_file)
}
aspect_file = file.path(morph_folder, "aspect.asc")
slope_file = file.path(morph_folder, "slope.asc")
fdir_file = file.path(morph_folder, "fdir.asc")
facc_file = file.path(morph_folder, "facc.asc")

dem = dem_file %>% rast %>% mask(roi_mask)
aspect = aspect_file %>% rast %>% mask(roi_mask)
slope = slope_file %>% rast %>% mask(roi_mask)
fdir = fdir_file %>% rast %>% mask(roi_mask)
facc = facc_file %>% rast %>% mask(roi_mask)


writeRaster(dem, file.path(morph_folder, "dem.asc"), overwrite = TRUE, NAflag = -9999)
writeRaster(aspect, file.path(morph_folder, "aspect.asc"), overwrite = TRUE, NAflag = -9999)
writeRaster(slope, file.path(morph_folder, "slope.asc"), overwrite = TRUE, NAflag = -9999)
writeRaster(fdir, file.path(morph_folder, "fdir.asc"), overwrite = TRUE, NAflag = -9999, datatype = "INT4S")
writeRaster(facc, file.path(morph_folder, "facc.asc"), overwrite = TRUE, NAflag = -9999, datatype = "INT4S")

