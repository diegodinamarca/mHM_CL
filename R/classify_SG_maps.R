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
                "gridExtra",
                "raster",
                "rasterVis",
                "ggpointdensity",
                "viridis",
                "scico") # list of packages to load
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


# Match with mhm resolution and extent ------------------------------------

outdir = c("DATA/SOIL_DATA/SoilGrids/soil_classification/")
dir.create(outdir)
r = rast("D:/mHM/DATA/SOIL_DATA/SoilGrids/SoilGrids_bd_clay_sand_0_200.tif")
plot(r)

# refernce morph data
reproject_to_mhm_morph <- function(r, morph.dir = "BASIN_DATA/data/7339001/morph_cr2soil",
                                   plot = FALSE
) {
  files = list.files(morph.dir, pattern = "asc$", full.names = TRUE)
  r2 = rast(files);r2
  if(plot){
    plot(r2)
  }
  project(r, r2, method = "bilinear")
}

r = reproject_to_mhm_morph(r, plot = FALSE)
plot(r[[1]])

# Rellenar pixeles NA
!is.na(r) %>% sum %>% plot


focal_repeat = function(r, n){
  if (n == 0){
    r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
    return(focal_repeat(r, n-1))
  }else{
    return(r)
    }
}
rf = focal_repeat(r, 3)
plot(rf[[1]])

dir.create(str_c(outdir, "NA_filled"))
writeRaster(rf, filename = paste0(outdir, "NA_filled/SoilGrids_bd_clay_sand_0_200.tif"))

# Create mhm classification maps ------------------------------------------
r = rast("D:/mHM/DATA/SOIL_DATA/SoilGrids/NA_filled/SoilGrids_bd_clay_sand_0_200.tif")
r %>% plot
depths = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")

# df to hold unique combinations from all the layers
data.depths <- data.frame()

# imagenes por profundidad
get_unique_combinations <- function(r, depths, var_names = c("da","arcilla","arena")) {
  depths = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")
  lapply(depths, function(d){
    cat("Cargando mapas de profundidad ", d,"\n")
    index = grep(d, names(r))
    imgs = r[[index]]
    cat("Extrayendo valores\n")
    imgs.df <- values(imgs)
    colnames(imgs.df) <- var_names
    imgs.df <- imgs.df %>% as_tibble %>% 
      mutate(arcilla = round(arcilla/10, digits = 0),
             arena = round(arena/10, digits = 0),
             da = round(da/100, digits = 1))
  }) %>% 
    bind_rows() %>% 
    group_by_all() %>%
    slice(1) %>%
    ungroup()

}
r
data.depths = get_unique_combinations(r)
glimpse(data.depths)

# data.depths %>% group_by_all() %>% slice(1)

# contar combinaciones presentes
comb.data <- data.depths %>% group_by_all() %>% count

# definir ID de tipo de suelo
comb.data$comb_id <- seq(1, nrow(comb.data));comb.data
# beepr::beep(sound = 5)

for (i in 1:6) {
  # i=1
  d <- depths[i]
  # cargar imagenes
  cat("Cargando mapas de profundidad ", d,"\n")
  index = grep(d, names(r))
  imgs = r[[index]]
  # extraer valores
  cat("Extrayendo valores...\n")
  imgs.df <- values(imgs)
  colnames(imgs.df) <- c("da","arcilla","arena")
  
  # redondear texturas a enteros y Da a 1 decimal
  imgs.df = imgs.df %>% 
    as_tibble() %>% 
    mutate(arcilla = round(arcilla/10, digits = 0),
           arena = round(arena/10, digits = 0),
           da = round(da/100, digits = 1))
  #asignar el ID de suelo a cada pixel
  cat("Asignando clases de suelo...\n")
  imgs.class <- left_join(imgs.df, comb.data, by = c("arcilla","arena","da"), keep = F)
  imgs.class
  #asignar valores al raster
  v <- rast(r[[1]])
  values(v) <- imgs.class$comb_id
  names(v) = paste0("soilclass_0",i)

  # png(paste0(outdir,"SG_class_horizon_0",i,".png"), width = 1400, height = 800)
  # plot(v)
  # dev.off()
  
  #exportar resultado
  cat("Exportando resultados...\n")
  names(v) = paste0("soil_class_horizon_0",i,".tif")
  v %>% writeRaster(str_c(outdir,"SG_class_horizon_0",i,".tif"), overwrite = T)
}



# exportar tabla de clasificaciones
comb.data %>% write_csv(paste0(outdir,"soilclass_withNA.csv"))

# ruta de imagenes
imgs <- list.files(outdir, full.names = T, pattern = "tif$") %>% rast
files <- list.files(outdir, full.names = T, pattern = "tif$");files

# Remover clases con valores NA
rm.class <- comb.data %>% ungroup() %>% filter(is.na(arena) | is.na(arcilla) | is.na(da)) %>%
  select(comb_id) %>% unique
rm.class <- rm.class$comb_id
rm.class

# exportar clases sin valores NA
comb.data %>% 
  filter(!(comb_id %in% rm.class)) %>%
  write_csv(paste0(outdir,"soilclasses.csv"))

EX.TABLE = comb.data %>% 
  mutate(ID = cur_group_id()) %>% 
  select(
    ID,
    `CLAY[%]` = arcilla,
    `SAND[%]` = arena,
    `Bd_mu[gcm-3]` = da,
    nSamples = comb_id
  )

nSoiltypes = nrow(comb.data);nSoiltypes
write_lines(str_c("nSoil_Types\t" ,as.integer(nSoiltypes)) ,str_c(outdir, "soil_classdefinition_iFlag_soilDB_1.txt"))
write_delim(EX.TABLE, str_c(outdir, "soil_classdefinition_iFlag_soilDB_1.txt"), 
            append = TRUE, 
            col_names = TRUE, delim = "\t")



# transform to ASCII
rr <- rast(imgs)
rr
names(imgs)
for (i in 1:6) {
  # i = 1
  print(i)
  img <- r[[i]]
  img[img %in% rm.class] <- NA
  rr <- c(rr, img)
}
is.na(rr) %>% plot
values(rr)

# beepr::beep(sound = 5)
dir.create(str_c(outdir, "NA_removed"))
names(rr) = paste0("soil_class_horizon_0",c(1,2,3,4,5,6),".tif")
writeRaster(rr, 
            paste0(outdir,"NA_removed/soil_class_horizon_0",c(1,2,3,4,5,6),".tif"),
            overwrite = T)
writeRaster(rr, 
            paste0(outdir,"NA_removed/soil_class_horizon_0",c(1,2,3,4,5,6),".asc"),
            overwrite = T, datatype = "INT4S")

