library(terra)
library(tidyverse)
library(stars)
library(cowplot)
source("R/utils.R")

theme_grid = theme(
  panel.background = element_rect(fill = "snow", ),
  panel.grid = element_line(color = "grey"),
  plot.background = element_rect(fill = "white"))


# Mosaicked outputs 
files.mosaic = list.files("../domain_chile/OUT", full.names = TRUE, pattern = ".nc$")
var.names = str_sub(basename(files.mosaic), end = -4)

# Domain folders to look for meteo forcings in monthly timesteps
domain_folders = dir("..", pattern = "domain_zone",  full.names = TRUE)
in.folder = file.path(domain_folders, "out/meteo")
files.domain = list.files(in.folder, full.names = TRUE) %>%
  list.files(full.names = TRUE, pattern = "annual")
files.domain


r.list = lapply(files.mosaic, function(f){
  r = rast(f)
  vn = varnames(r)
  if (vn %in% c("aET","Q")){
    r_ann = annual_mean(r, "sum")
  }else{
    r_ann = annual_mean(r, "mean")
  }
})
names(r.list) = var.names

# load and mosaic meteo forcings
pet.files = files[grep("pet", files.domain)]
pre.files = files[grep("pre", files.domain)]

pet.list = lapply(pet.files, rast)
pet = do.call(mosaic, pet.list)
varnames(pet) = "pet"

pre.list = lapply(pre.files, rast)
pre = do.call(mosaic, pre.list)
varnames(pre) = "pre"

# add to list of rasters
r.list$pre = pre
r.list$pet = pet
r.list$`P-ET` = r.list$pre - r.list$aET
r.list
length(r.list)

# Order and variable info
order_top <- c("pre", "pet", "aET", "P-ET")
order_bottom <- c("Q", "SM_Lall", "satSTW", "snowpack")

# Helper function to plot one raster from the list
plot_from_list <- function(varname) {
  r <- r.list[[varname]]
  if (varname == "SM_Lall"){
    limits <- get_scale_lim(r, round.digit = 2)
  }else{
    limits <- get_scale_lim(r)
  }
  plot_raster(r, depth = "", var = str_c(varname, "[mm]"), limits = limits)
}

# Generate plots
plots_top <- lapply(order_top, plot_from_list)
plots_bottom <- lapply(order_bottom, plot_from_list)

# Arrange the plots in a grid
final_plot <- ggpubr::ggarrange(
  plotlist = c(plots_top, plots_bottom),
  ncol = 4, nrow = 2,
  common.legend = FALSE
)

# Display the plot
print(final_plot+theme(plot.background = element_rect(fill = "white")))
ggsave(filename = file.path("../FIGS","annual_mean_output_forcings.png"), width = 7)
       