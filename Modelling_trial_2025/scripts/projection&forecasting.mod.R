library(terra)
library(ggplot2)
library(tidyterra)


## ---- CURRENT PROJECTIONS ----
# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  #metric.binary = 'all',
                                  #metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)


# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'CurrentEM',
                                             new.env = myExpl,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)


## ---- WHAT ABOUT THE FUTURE PROJECTIONS? -----

# Prepare Bio-Oracle layers for future (2020-2100)
layers.bio2 <- list_layers(simplify=T) #simplify=F if need more info

# SSP1-1.9 - Very low emissions scenario, aligned with the Paris Agreement goal of limiting global warming to 1.5°C above pre-industrial levels.
# SSP1-2.6 - A low-emissions scenario, limiting warming to around 2°C by 2100.
# SSP2-4.5 - Comparable to current policy trends. A "middle-of-the-road" scenario, with moderate emissions and warming around 2.5–3°C by 2100.
# SSP3-7.0 - A high-emissions scenario, leading to around 4°C warming by 2100.
# SSP4-6.0 - A moderate-to-high emissions scenario, with warming of around 3°C by 2100.
# SSP5-8.5 - Very high-emissions scenario, leading to more than 4°C warming by 2100.

# Load required libraries
library(terra)

# Define the base directory for SSP layers
#base_dir <- "C:/Users/Justine/OneDrive/Documents/PhD/biomod2_git/layers/future"
base_dir <- "/home/gunnarandersson/Dokument/Projekt/HAV2025/data/Biooracle.download/future"

if (!dir.exists(base_dir)) dir.create(base_dir) 

# Define SSP scenarios and their corresponding dataset lists
ssp_scenarios <- list(
  "SSP1" = list(
    datasets = list(
      list(dataset_id = "o2_SSP1_2020_2100_depthmean", variables = c("o2_mean")),
      list(dataset_id = "chl_SSP1_2020_2100_depthmean", variables = c("chl_mean")),
      list(dataset_id = "thetao_SSP1_2020_2100_depthmean", variables = c("thetao_mean")),
      list(dataset_id = "ph_SSP1_2020_2100_depthmean", variables = c("ph_mean")),
      list(dataset_id = "po4_SSP1_2020_2100_depthmean", variables = c("po4_mean")),
      list(dataset_id = "so_SSP1_2020_2100_depthmean", variables = c("so_mean")),
      list(dataset_id = "phyc_SSP1_2020_2100_depthmean", variables = c("phyc_mean"))
    )
  ),
  
  "SSP2" = list(
    datasets = list(
      list(dataset_id = "o2_SSP2_2020_2100_depthmean", variables = c("o2_mean")),
      list(dataset_id = "chl_SSP2_2020_2100_depthmean", variables = c("chl_mean")),
      list(dataset_id = "thetao_SSP2_2020_2100_depthmean", variables = c("thetao_mean")),
      list(dataset_id = "ph_SSP2_2020_2100_depthmean", variables = c("ph_mean")),
      list(dataset_id = "po4_SSP2_2020_2100_depthmean", variables = c("po4_mean")),
      list(dataset_id = "so_SSP2_2020_2100_depthmean", variables = c("so_mean")),
      list(dataset_id = "phyc_SSP2_2020_2100_depthmean", variables = c("phyc_mean"))
    )
  ),
  
  "SSP3" = list(
    datasets = list(
      list(dataset_id = "o2_SSP3_2020_2100_depthmean", variables = c("o2_mean")),
      list(dataset_id = "chl_SSP3_2020_2100_depthmean", variables = c("chl_mean")),
      list(dataset_id = "thetao_SSP3_2020_2100_depthmean", variables = c("thetao_mean")),
      list(dataset_id = "ph_SSP3_2020_2100_depthmean", variables = c("ph_mean")),
      list(dataset_id = "po4_SSP3_2020_2100_depthmean", variables = c("po4_mean")),
      list(dataset_id = "so_SSP3_2020_2100_depthmean", variables = c("so_mean")),
      list(dataset_id = "phyc_SSP3_2020_2100_depthmean", variables = c("phyc_mean"))
    )
  ),
  
  "SSP4" = list(
    datasets = list(
      list(dataset_id = "o2_SSP4_2020_2100_depthmean", variables = c("o2_mean")),
      list(dataset_id = "chl_SSP4_2020_2100_depthmean", variables = c("chl_mean")),
      list(dataset_id = "thetao_SSP4_2020_2100_depthmean", variables = c("thetao_mean")),
      list(dataset_id = "ph_SSP4_2020_2100_depthmean", variables = c("ph_mean")),
      list(dataset_id = "po4_SSP4_2020_2100_depthmean", variables = c("po4_mean")),
      list(dataset_id = "so_SSP4_2020_2100_depthmean", variables = c("so_mean")),
      list(dataset_id = "phyc_SSP4_2020_2100_depthmean", variables = c("phyc_mean"))
    )
  ),
  
  "SSP5" = list(
    datasets = list(
      list(dataset_id = "o2_SSP5_2020_2100_depthmean", variables = c("o2_mean")),
      list(dataset_id = "chl_SSP5_2020_2100_depthmean", variables = c("chl_mean")),
      list(dataset_id = "thetao_SSP5_2020_2100_depthmean", variables = c("thetao_mean")),
      list(dataset_id = "ph_SSP5_2020_2100_depthmean", variables = c("ph_mean")),
      list(dataset_id = "po4_SSP5_2020_2100_depthmean", variables = c("po4_mean")),
      list(dataset_id = "so_SSP5_2020_2100_depthmean", variables = c("so_mean")),
      list(dataset_id = "phyc_SSP5_2020_2100_depthmean", variables = c("phyc_mean"))
    )
  )
)

# Ensure directories exist for each SSP scenario
ssp_dirs <- file.path(base_dir, names(ssp_scenarios))
lapply(ssp_dirs, function(dir) if (!dir.exists(dir)) dir.create(dir))

# Define the time constraint for all scenarios
constraints <- list(
  time = c("2020-01-01T00:00:00Z", "2100-01-01T00:00:00Z")  # Adjust time constraints as needed
)

# Loop through each SSP scenario and process the layers
ssp_layers <- list()

for (ssp in names(ssp_scenarios)) {
  
  message("Processing ", ssp)
  
  ssp_dir <- file.path(base_dir, ssp) # Define the directory for the current scenario
  datasets <- ssp_scenarios[[ssp]]$datasets # Retrieve dataset list for the SSP
  
  # Download datasets
  for (dataset in datasets) {
    download_layers(dataset$dataset_id, variables = dataset$variables, constraints = constraints, directory = ssp_dir)
  }
  
  # Load all downloaded raster layers for this SSP
  layers <- rast(paste0(ssp_dir, "/", list.files(ssp_dir)))
  
  # Store layers in a list for future use
  ssp_layers[[ssp]] <- layers
  
  # Plot layers
  for (i in 1:nlyr(layers)) {
    plot(layers[[i]], main = paste(ssp, "-", names(layers)[i]))
  }
  
  # Assign layers to biomod2 variable dynamically
  assign(paste0("myExplFuture", ssp), layers, envir = .GlobalEnv)
}

# Print summary of loaded layers
print(ssp_layers)


# Project onto future conditions
myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                              proj.name = 'Future',
                                              new.env = myExplFuture,
                                              models.chosen = 'all',
                                              metric.binary = 'TSS',
                                              build.clamping.mask = TRUE)


## ---- Compare range size for current and future predictions ----
# Load current and future binary projections
CurrentProj <- get_predictions(myBiomodProj, metric.binary = "TSS")
FutureProj <- get_predictions(myBiomodProjectionFuture, metric.binary = "TSS")

# Compute differences
myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = CurrentProj, 
                                      proj.future = FutureProj)

myBiomodRangeSize$Compt.By.Models
plot(myBiomodRangeSize$Diff.By.Pixel)

# Represent main results 
gg = bm_PlotRangeSize(bm.range = myBiomodRangeSize, 
                      do.count = TRUE,
                      do.perc = TRUE,
                      do.maps = TRUE,
                      do.mean = TRUE,
                      do.plot = TRUE,
                      row.names = c("Species", "Dataset", "Run", "Algo"))


