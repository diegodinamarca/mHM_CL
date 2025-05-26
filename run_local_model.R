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
  # Script Description:
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
                "scico", 
                "shiny", 
                "rlist") # list of packages to load
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
  basin_mean <- function(flux) {
    as_tibble(values(flux)) %>% 
      pivot_longer(cols = 1:ncol(.), names_to = "time", values_to = "var") %>% 
      mutate(time = ymd(time)) %>%
      drop_na() %>%
      group_by(time) %>%
      summarise(var = mean(var))
  }
  
  basin_mean_bylayer <- function(flux) {
    data_day = lapply(1:6, function(i){
      if (!is.null(flux[[i]])){
        x = as_tibble(values(flux[[i]])) %>% 
          pivot_longer(cols = 1:nlyr(flux[[i]]), names_to = "time",values_to = "var") %>% 
          drop_na() %>% 
          mutate(time = ymd(time)) %>% 
          group_by(time) %>%
          summarise(var = mean(var))
      }
    }) %>% bind_cols()
    data_day = data_day %>% 
      select(time = 1, starts_with("var"))
  }
}


models = c("default.2m","SG.calib.2m", "FAO.calib.2m","default.4m","SG.calib.4m", "FAO.calib.4m")
prev.results = TRUE
climdir =  "meteo_cr2met_v2.5_1960_2023"
input.file = here("PROC","mHM_basin_mean_day.csv")
output.file = here("PROC","mHM_basin_mean_day_update.csv")
default.model = "default.2m"

# input$outputs = names(output.list)

# Define the outputs
output.list= c("Total aET of the soil profile",
               "Actual ET by soil layer",
               "Potential evapotranspiration (mHM corrected)",
               "Infiltration",
               "Percolation",
               "Soil Moisture [v/v]",
               "Soil Water content [mm]",
               "Total Runoff",
               "Rapid Interflow",
               "Slow Interflow",
               "Base flow")

names(output.list) = c("ETR","aET","PET","inf","per","sm","swc","roff","rflow","sflow","bflow")
# input = list()
# input$basin_id = 7339001
# input$models_to_process = models
model_results = list()

# Define the user interface
ui <- fluidPage(
  
  # App title
  titlePanel("Model Processing App"),
  
  # Sidebar layout for inputs
  sidebarLayout(
    sidebarPanel(
      
      # Numeric input for basin id
      numericInput("basin_id", "Enter Basin ID:", value = 7339001, min = 1),  # Set default value as needed
      
      # Checkbox input for process.forcing
      checkboxInput("process_forcing", "Process Forcings of the models", value = TRUE),

      # Dynamic model selection, based on models.update
 
      checkboxGroupInput("models_to_process", "Select Models to Process:",
                         choices = models),   # Assuming 'models' is a vector in your script
      # Section to select outputs to be processed
      checkboxGroupInput("outputs", "Select Outputs to Process:",
                         choices = setNames(names(output.list), output.list)),
      
      # Action button to start processing
      actionButton("process", "Run Processing")
    ),
    
    # Main panel for displaying messages
    mainPanel(
      verbatimTextOutput("console_output")  # This will display processing messages
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  
  # Store messages
  messages <- reactiveVal("")
  
  # Append message function
  append_message <- function(msg) {
    current_messages <- messages()
    messages(paste(current_messages, msg, sep = "\n"))
  }
  
  # Observe when the process button is clicked
  observeEvent(input$process, {
    
    # Clear previous messages
    messages("")
    
    # Display the entered basin ID
    append_message(paste("Processing for Basin ID:", input$basin_id))
    
    append_message("Loading helper functions")
    source(here::here("scripts","helper_functions.R"))
    
    append_message("Loading basin limits")
    cuenca = read_sf("DATA/SHP/Cauquenes.shp") %>% st_transform(4326)
    
    
    # # Checking for previous output file
    # if (file.exists(input.file)){
    #   data_day = read_csv(input.file) %>%
    #     select(-starts_with("q"))
    #   data_day %>% names
    #   
    #   for (i in 1:length(input$models_to_process)) {
    #     
    #       col.exists = str_match(names(data_day), input$models_to_process[i]) %>%
    #         as.factor %>%
    #         as.numeric %>%
    #         sum(na.rm=TRUE)
    #       
    #       if (col.exists){
    #         data_day = data_day %>% select(-ends_with(input$models_to_process[i]))
    #         message("removed columns: ",models[i])
    #       }
    #     
    #   }
    #   process.default = TRUE
    #   if (process.default){
    #     col.exists = str_match(names(data_day), "def") %>%
    #       as.factor %>%
    #       as.numeric %>%
    #       sum(na.rm=TRUE)
    # 
    #     if (col.exists){
    #       data_day = data_day %>% select(-ends_with("def"))
    #       message("removed columns: ","def")
    #     }
    #   }
    # }
    # else {
    #   prev.results = FALSE
    # }
    # 
    # Start processing models
    n = length(input$models_to_process)
    if (n != 0) {
      for (model in input$models_to_process) {
        append_message(paste("Processing model:", model))

        # process.default = TRUE
        if (!is.null(input$outputs)) {
          if ("ETR" %in% input$outputs) {
            append_message("Processing ETR output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 10, daily = TRUE, process = FALSE)
            model_results[[model]]$ET = basin_mean(flux)
            names(model_results[[model]]$ET) = c("time", str_c("ETR_",model))
            # Insert code for processing fluxes here
          }
          if ("PET" %in% input$outputs) {
            append_message("Processing PET output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 9, daily = TRUE, process = FALSE)
            model_results[[model]]$PET = basin_mean(flux)
            names(model_results[[model]]$PET) = c("time", str_c("PET_",model))
          }
          if ("per" %in% input$outputs) {
            append_message("Processing percolation output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 16, daily = TRUE, process = FALSE)
            model_results[[model]]$per = basin_mean(flux)
            names(model_results[[model]]$per) = c("time", str_c("per_",model))
          }
          if ("roff" %in% input$outputs) {
            append_message("Processing total runoff output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 11, daily = TRUE, process = FALSE)
            model_results[[model]]$roff = basin_mean(flux)
            names(model_results[[model]]$roff) = c("time", str_c("roff_",model))
          }
          if ("rflow" %in% input$outputs) {
            append_message("Processing rapid flow output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 13, daily = TRUE, process = FALSE)
            model_results[[model]]$rflow = basin_mean(flux)
            names(model_results[[model]]$rflow) = c("time", str_c("rflow_",model))
          }
          if ("sflow" %in% input$outputs) {
            append_message("Processing slow flow output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 14, daily = TRUE, process = FALSE)
            model_results[[model]]$sflow = basin_mean(flux)
            names(model_results[[model]]$sflow) = c("time", str_c("sflow_",model))
          }
          if ("bflow" %in% input$outputs) {
            append_message("Processing base flow output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 14, daily = TRUE, process = FALSE)
            model_results[[model]]$bflow = basin_mean(flux)
            names(model_results[[model]]$bflow) = c("time", str_c("bflow_",model))
          }
          if ("aET" %in% input$outputs) {
            append_message("Processing ET by layer output.")
            # Insert code for processing water balance here
            flux = get_flux(model, idbasin = input$basin_id, n = 19, daily = TRUE, process = FALSE)
            df.m = basin_mean_bylayer(flux)
            names(df.m)[2:ncol(df.m)] = paste0("aET_L", 1:(ncol(df.m)-1),"_",model)
            model_results[[model]]$aET = df.m
            names(model_results[[model]]$aET)[2:7] = paste0("aET_L", 1:6,"_",model)
          }
          if ("sm" %in% input$outputs) {
            append_message("Processing soil moisture [v/v] output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 4, daily = TRUE, process = FALSE)
            df.m = basin_mean_bylayer(flux)
            names(df.m)[2:ncol(df.m)] = paste0("sm_L", 1:(ncol(df.m)-1),"_",model)
            model_results[[model]]$sm = df.m
            names(model_results[[model]]$sm)[2:7] = paste0("sm_L", 1:6,"_",model)
          }
          if ("swc" %in% input$outputs) {
            append_message("Processing soil moisture [mm] output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 3, daily = TRUE, process = FALSE)
            df.m = basin_mean_bylayer(flux)
            names(df.m)[2:ncol(df.m)] = paste0("swc_L", 1:(ncol(df.m)-1),"_",model)
            model_results[[model]]$swc = df.m
            names(model_results[[model]]$swc)[2:7] = paste0("swc_L", 1:6,"_",model)
          }
          if ("inf" %in% input$outputs) {
            append_message("Processing infiltration output.")
            flux = get_flux(model, idbasin = input$basin_id, n = 17, daily = TRUE, process = FALSE)
            df.m = basin_mean_bylayer(flux)
            names(df.m)[2:ncol(df.m)] = paste0("inf_L", 1:(ncol(df.m)-1),"_",model)
            model_results[[model]]$aET = df.m
            names(model_results[[model]]$aET)[2:7] = paste0("inf_L", 1:6,"_",model)
          }
        } else {
          append_message("No outputs selected for processing.")
        }
        
      } 
    } else{
      append_message("Error. Please select the models you want to process")
    }
    
    # Process based on forcing/default options
    if (input$process_forcing) {
      append_message("Processing with forcing enabled.")
      # process.forcing = TRUE
      pp = rast(here("BASIN_DATA","DATA",input$basin_id,climdir,"pre.NC"))
      names(pp) = time(pp)
      pp_day = extract(pp, cuenca, fun = "mean") %>% 
        as_tibble() %>% 
        pivot_longer(cols = 2:ncol(.), names_to = "time",values_to = "PP") %>% 
        mutate(time = ymd(time)) %>% 
        select(-ID)
      model_results$clim = pp_day
      
      #   PET ---------------------------------------------------------------------
      dates <- seq(ymd("1960-01-01"),ymd("2021-12-31"), by = "days");length(dates)
      pet = rast(here("BASIN_DATA","DATA",input$basin_id,climdir,"pet.NC"))
      names(pet) <- time(pet)
      pet_day = extract(pet, cuenca, fun = "mean") %>% 
        as_tibble() %>% 
        pivot_longer(cols = 2:ncol(.), names_to = "time",values_to = "PET_cr2") %>% 
        mutate(time = ymd(time)) %>% 
        select(-ID)
      model_results$clim = full_join(pp_day, pet_day, by = "time")
      
    } else {
      append_message("Skipping forcing processing.")
      # process.forcing = FALSE
    }
    
    append_message("Processing complete.")
    # time = output[1]$time
    
    write_rds(model_results, here(str_c(input$basin_id, "_model_results.rds")))
    append_message(str_c("Results saved to ",here(str_c(input$basin_id, "_model_results.rds"))))
    
    # Output console messages
    output$console_output <- renderText({
      messages()
    })
  } )
  
}

# Run the app
shinyApp(ui = ui, server = server)
