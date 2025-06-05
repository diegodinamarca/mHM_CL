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
  source("scripts/R/utils.R")
}


# Match with mhm resolution and extent ------------------------------------
config <- fromJSON("new_domain/preprocess_config.json")

roi_file = config$roi_file
soil_files = config$soil_folder
cellsize = config$cellsize
ref_file = file.path(domain_path, config$meteo_folder,"pre.nc")
morph_folder = file.path(domain_path, config$morph_folder)
temp_folder = file.path(morph_folder, "temp");temp_folder
dir.create(temp_folder, showWarnings = FALSE)

# read ROI
roi = read_sf(roi_file)
box = get_extent(roi)
# Crear buffer de 1km (3000 metros)
roi_buff <- st_buffer(box, dist = 3000)
# read reference grid
ref = rast(ref_file)
res(ref) = cellsize
# refernce morph data
reproject_to_mhm_morph <- function(r) {
  r = crop(r, roi_buff)
  # r = crop(r, ref)
  r = resample(r, ref, method = "bilinear")
}

files = list.files(soil_files, full.names = TRUE)
bulkd_files = files[grep("Bulkd", files)]
clay_files = files[grep("Clay", files)]
sand_files = files[grep("Sand", files)]

bulkd_r = reproject_to_mhm_morph(rast(bulkd_files))
clay_r = reproject_to_mhm_morph(rast(clay_files))
sand_r = reproject_to_mhm_morph(rast(sand_files))


focal_repeat = function(r, n){
  if (n != 0){
    r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
    return(focal_repeat(r, n-1))
  }else{
    return(r)
    }
}

bulkd_filled = focal_repeat(bulkd_r, 10)
plot(bulkd_r[[1]])
plot(bulkd_filled[[1]])
writeRaster(bulkd_filled, filename = paste0(temp_folder, "/bulkd_filled.tif"), overwrite = TRUE)

clay_filled = focal_repeat(clay_r, 10)
plot(clay_r[[1]])
plot(clay_filled[[1]])
writeRaster(clay_filled, filename = paste0(temp_folder, "/clay_filled.tif"), overwrite = TRUE)

sand_filled = focal_repeat(sand_r, 10)
plot(sand_r[[1]])
plot(sand_filled[[1]])
writeRaster(sand_filled, filename = paste0(temp_folder, "/sand_filled.tif"), overwrite = TRUE)

# Create mhm classification maps ------------------------------------------
depths = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")

# df to hold unique combinations from all the layers
data.depths <- data.frame()

# imagenes por profundidad
get_combinations <- function(r, var_names = c("Bulkd","Clay","Sand")) {
  lapply(depths, function(d){
    cat("Cargando mapas de profundidad ", d,"\n")
    index = grep(d, names(r))
    imgs = r[[index]]
    cat("Extrayendo valores\n")
    imgs.df <- values(imgs)
    colnames(imgs.df) <- var_names
    imgs.df %>% as_tibble()
  }) %>% 
    bind_rows()
}

all_soil = c(
  round(bulkd_filled,digits=2),
  round(sand_filled,digits=0),
  round(clay_filled,digits=0)
  )

soil.mask = sum(!is.na(all_soil))
soil.mask[soil.mask != 18] = NA
all_soil = mask(all_soil, soil.mask)

data.depths = get_combinations(all_soil)
glimpse(data.depths)

# contar combinaciones presentes
comb.data <- data.depths %>% group_by_all() %>% count() %>% drop_na()

# definir ID de tipo de suelo
comb.data$comb_id <- seq(1, nrow(comb.data));comb.data

# beepr::beep(sound = 5)

for (i in 1:6) {
  # i=1
  d <- depths[i]
  # cargar imagenes
  cat("Cargando mapas de profundidad ", d,"\n")
  index = grep(d, names(all_soil))
  imgs = all_soil[[index]]
  # extraer valores
  cat("Extrayendo valores...\n")
  imgs.df <- values(imgs)
  colnames(imgs.df) <- c("Bulkd","Clay","Sand")
  
  # redondear texturas a enteros y Da a 1 decimal
  imgs.df = imgs.df %>% 
    as_tibble()
  #asignar el ID de suelo a cada pixel
  cat("Asignando clases de suelo...\n")
  imgs.class <- left_join(imgs.df, comb.data, by = c("Clay","Sand","Bulkd"), keep = F)
  imgs.class
  #asignar valores al raster
  v <- rast(all_soil[[1]])
  values(v) <- imgs.class$comb_id
  names(v) = paste0("soilclass_0",i)

  #exportar resultado
  cat("Exportando resultados...\n")
  names(v) = paste0("soil_class_horizon_0",i,".tif")
  v %>% writeRaster(str_c(temp_folder,"/soil_class_horizon_0",i,".tif"), overwrite = T)
}

# ruta de imagenes
imgs <- list.files(temp_folder, full.names = T, pattern = "soil_class_horizon_[0-9]+.tif$") %>% rast;imgs

EX.TABLE = comb.data %>% 
  mutate(ID = comb_id) %>% 
  mutate(Clay = as.integer(Clay),
         Sand = as.integer(Sand)) %>% 
  select(
    ID,
    `CLAY[%]` = Clay,
    `SAND[%]` = Sand,
    `Bd_mu[gcm-3]` = Bulkd,
    nSamples = n
  )


nSoiltypes = nrow(comb.data);nSoiltypes
write_lines(str_c("nSoil_Types\t" ,as.integer(nSoiltypes)) ,str_c(morph_folder, "/soil_classdefinition_iFlag_soilDB_1.txt"))
write_delim(EX.TABLE, str_c(morph_folder, "/soil_classdefinition_iFlag_soilDB_1.txt"), 
            append = TRUE, 
            col_names = TRUE, delim = "\t")

# transform to ASCII
rr <- rast(imgs)
names(imgs)
for (i in 1:6) {
  # i = 1
  print(i)
  img <- imgs[[i]]
  rr <- c(rr, img)
}
is.na(rr) %>% plot

# beepr::beep(sound = 5)
names(rr) = paste0("soil_class_horizon_0",c(1,2,3,4,5,6),".tif")
writeRaster(rr,
            paste0(morph_folder,"/soil_class_horizon_0",c(1,2,3,4,5,6),".asc"),
            overwrite = T, datatype = "INT4S", NAflag = -9999)

