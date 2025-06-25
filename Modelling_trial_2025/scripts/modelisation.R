#install.packages("biomod2")
library(biomod2)
library(terra)
library(raster)
library(devtools)
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(rgbif)
library(sf)
library(dplyr)
library(ggplot2)
library(maps)

setwd("C:/Users/Justine/OneDrive/Documents/PhD/biomod2_git") 

# --------- Data preparation ----------

## OUR SPECIES OCCURRENCES ---
# Load species occurrences
# Columns: Lat, Lon, Species Names...
# Rows: coordinates, P/A data (0 and 1)

DataSpecies <- read.csv2("data_test.csv",header=TRUE, sep=",",stringsAsFactors = F)
DataSpecies_bySite <- DataSpecies %>%
  group_by(Site) %>%
  summarise(
    longitude = first(longitude),  # Keep first longitude value per site
    latitude = first(latitude),    # Keep first latitude value per site
    across(where(is.numeric) & !c(longitude, latitude), max, na.rm = TRUE) # Apply max to species columns
  )
write.table(DataSpecies,"ARMS_NIS_bySite.csv")

DataSpecies_byObs <- DataSpecies %>%
  group_by(Observatory) %>%
  summarise(
    longitude = first(longitude),  # Keep first longitude value per site
    latitude = first(latitude),    # Keep first latitude value per site
    across(where(is.numeric) & !c(longitude, latitude), max, na.rm = TRUE) # Apply max to species columns
  )
write.table(DataSpecies_byObs,"ARMS_NIS_byObs.csv")

# Select the name of the studied species
myRespName <- 'Fenestrulinadelicia'

# Get corresponding presence/absence data
myResp <- DataSpecies_byObs[, myRespName]
myResp <- as.numeric(unlist(myResp))
table(myResp)  # Check presence-absence distribution

# Get lat/lon 
# CAREFUL! Some points might be too close to the shore to obtain env data
myRespXY <- DataSpecies_byObs[, c('longitude', 'latitude')]
myRespXY <- as.data.frame(myRespXY)
myRespXY$longitude <- as.numeric(myRespXY$longitude)
myRespXY$latitude <- as.numeric(myRespXY$latitude)

# Vizualize the points to be sure
world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = myRespXY, aes(x = longitude, y = latitude), color = "red", alpha = 0.7) +
  labs(title = "Sampling Locations", x = "Longitude", y = "Latitude") +
  coord_cartesian(xlim = c(-15,40), ylim = c(25,80))
  theme_minimal()


## BIO-ORACLE LAYERS
# Explore environmental variables from Bio-Oracle
layers.bio2 <- list_layers(simplify=F) #simplify=F if need more info

# Download and load data layers of choice 
dir <- "C:/Users/Justine/OneDrive/Documents/PhD/biomod2_git/layers"
if (!dir.exists(dir)) dir.create(dir) 

constraints <- list(
  #latitude = c(25, 80),
  #longitude = c(-15, 40),
  time = c("2010-01-01T00:00:00Z", "2010-01-01T00:00:00Z") # Ensure time is correctly specified
)

datasets <- list(
  
  list(dataset_id = "o2_baseline_2000_2018_depthmean",
       variables = c("o2_mean"),
       constraints = constraints),
  
  list(dataset_id = "chl_baseline_2000_2018_depthmean",
       variables = "chl_mean",
       constraints = constraints),
  
  list(dataset_id = "thetao_baseline_2000_2019_depthmean",
       variables = c("thetao_mean"),
       constraints = constraints),

  list(dataset_id = "ph_baseline_2000_2018_depthmean",
       variables = c("ph_mean"),
      constraints = constraints),
  
  list(dataset_id = "po4_baseline_2000_2018_depthmean",
       variables = c("po4_mean"),
       constraints = constraints),
  
  list(dataset_id = "so_baseline_2000_2019_depthmean",
      variables = c("so_mean"),
      constraints = constraints),
  
  list(dataset_id = "phyc_baseline_2000_2020_depthmean",
       variables = c("phyc_mean"),
       constraints = constraints))

datasets_depthmax <- list(
  
  list(dataset_id = "o2_baseline_2000_2018_depthmax",
     variables = c("o2_mean"),
     constraints = constraints),

  list(dataset_id = "chl_baseline_2000_2018_depthmax",
     variables = c("chl_mean"),
     constraints = constraints),

  list(dataset_id = "thetao_baseline_2000_2019_depthmax",
     variables = c("thetao_mean"),
     constraints = constraints),

  list(dataset_id = "ph_baseline_2000_2018_depthmax",
     variables = c("ph_mean"),
     constraints = constraints),

  list(dataset_id = "po4_baseline_2000_2018_depthmax",
     variables = c("po4_mean"),
     constraints = constraints),

  list(dataset_id = "so_baseline_2000_2019_depthmax",
     variables = c("so_mean"),
     constraints = constraints),

  list(dataset_id = "phyc_baseline_2000_2020_depthmax",
     variables = c("phyc_mean"),
     constraints = constraints))

datasets_depthmin <- list(
  
  list(dataset_id = "o2_baseline_2000_2018_depthmin",
       variables = c("o2_mean"),
       constraints = constraints),
  
  list(dataset_id = "chl_baseline_2000_2018_depthmin",
       variables = c("chl_mean"),
       constraints = constraints),
  
  list(dataset_id = "thetao_baseline_2000_2019_depthmin",
       variables = c("thetao_mean"),
       constraints = constraints),
  
  list(dataset_id = "ph_baseline_2000_2018_depthmin",
       variables = c("ph_mean"),
       constraints = constraints),
  
  list(dataset_id = "po4_baseline_2000_2018_depthmin",
       variables = c("po4_mean"),
       constraints = constraints),
  
  list(dataset_id = "so_baseline_2000_2019_depthmin",
       variables = c("so_mean"),
       constraints = constraints),
  
  list(dataset_id = "phyc_baseline_2000_2020_depthmin",
       variables = c("phyc_mean"),
       constraints = constraints))

for (dataset in c(datasets, datasets_depthmax)) {
  
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  download_layers(dataset_id, variables = variables, constraints = constraints, directory= dir)
}

layers <- rast(paste0(dir,"/", list.files(dir)[]))

for (i in 1:nlyr(layers)) {  # nlyr(layer) returns the number of layers in the raster
  plot(layers[[i]], main = names(layers)[i])  # Adds a title with the layer name
}

#Tentative merge with biomod2
myExpl <- layers # This is a SpatRaster
plot(myExpl)
names(myExpl)


# Define Europe bounding box (xmin, xmax, ymin, ymax)
europe_extent <- ext(-15, 40, 25, 80)

# Crop the raster to Europe
myExpl_europe <- crop(myExpl, europe_extent)

# Plot to verify
plot(myExpl_europe[[1]], main = "First Environmental Layer - Europe")




# little test of correlation between predictors, using the variance inflation factor (VIF)
library(usdm)
vifstep(layers)  # Check collinearity for all predictors
library(ENMTools)

# Pearson's correlation between layers. If a pair has a correlation factor >0.8, keep only the layer with the highest contribution
raster.cor.matrix(myExpl, method = "pearson")
raster.cor.plot(myExpl, method = "pearson")

#layers <- dropLayer(layers, which(names(layers) == ""))


### HERE WE SHOULD GET THE DATA FROM GBIF TO COMPLEMENT
species_name <- "Bonnemaisonia hamifera"  # Replace with your species name

# Query GBIF
occ_data <- occ_search(
  scientificName = species_name,
  hasCoordinate = TRUE,  # Only return geo-referenced occurrences
  limit = 10000)  # Increase if needed
  #geometry = "POLYGON((-15 25, 40 25, 40 80, -15 80, -15 25))" #only get europe otherwise it takes forever


gbif_occ <- occ_data$data
unique(gbif_occ$basisOfRecord)
gbif_occ <- gbif_occ[gbif_occ$basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE"), ]

gbif_occ_clean <- gbif_occ[, c("decimalLongitude", "decimalLatitude", "eventDate")]
colnames(gbif_occ_clean) <- c("longitude", "latitude", "date")
gbif_occ_clean$date <- as.Date(gbif_occ_clean$date)
gbif_occ_filtered <- gbif_occ_clean[
  gbif_occ_clean$date >= as.Date("2018-01-01") & gbif_occ_clean$date <= as.Date("2024-01-01"),
]
gbif_occ_filtered <- gbif_occ_filtered[, c("longitude", "latitude")]
gbif_occ_filtered <- unique(gbif_occ_filtered)
gbif_occ_filtered$presence <- 1 # Add GBIF occurrences as presence data


# OBIS data
library(robis)

#obtain occurrence data
moll<-occurrence("Mollusca")
moll_abs<-occurrence("Mollusca", absence="include") #include absence records
write.csv(moll, "mollusca-obis.csv")

#obtain a list of datasets
molldata<-dataset(scientificname="Mollusca")

#obtain a checklist of Mollusc species in a certain area
mollcheck<-checklist(scientificname="Mollusca", geometry = "POLYGON ((2.3 51.8, 2.3 51.6, 2.6 51.6, 2.6 51.8, 2.3 51.8))")



# Here we could set aside around 20% of GBIF data (randomly picked) to have a data set unseen by the model and use it as evaluation
# https://github.com/biomodhub/biomod2/issues/526


# Prepare my data before merge with GBIF data
my_data_clean <- data.frame(
  longitude = myRespXY$longitude,
  latitude = myRespXY$latitude,
  presence = myResp  # 1 for presence, 0 for absence
)

# Merge my data with GBIF
combined_data <- rbind(my_data_clean, gbif_occ_filtered)
combined_data <- unique(combined_data)
combined_data$presence <- as.numeric(combined_data$presence)
combined_data$longitude <- as.numeric(combined_data$longitude)
combined_data$latitude <- as.numeric(combined_data$latitude)

# Vizualize the points to be sure
world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = combined_data, aes(x = longitude, y = latitude), color = "red", alpha = 0.7) +
  labs(title = "Occurrences Locations", x = "Longitude", y = "Latitude") +
  #coord_cartesian(xlim = c(-15,40), ylim = c(25,80)) +
  theme_minimal()


world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = combined_data, 
             aes(x = longitude, y = latitude, color = factor(presence)), 
             alpha = 0.7) +
  scale_color_manual(values = c("0" = "red", "1" = "#789630"), 
                     labels = c("Absent", "Present"), 
                     name = "Presence") +
  labs(title = "Occurrences Locations", x = "Longitude", y = "Latitude") +
  #coord_cartesian(xlim = c(-15, 40), ylim = c(25, 80)) +
  theme_minimal()


# Get corresponding presence/absence data
myResp <- combined_data$presence
myResp <- as.numeric(unlist(myResp))
table(myResp)  # Check presence-absence distribution

# Get lat/lon 
# CAREFUL! Some points might be too close to the shore to obtain env data
myRespXY <- combined_data[, c('longitude', 'latitude')]
myRespXY <- as.data.frame(myRespXY)
myRespXY$longitude <- as.numeric(myRespXY$longitude)
myRespXY$latitude <- as.numeric(myRespXY$latitude)

# Here, go to shift_coordinates.R to correct the lat/lon of sampling points outside of env layers extent
# then put the corrected ones back in this script

myRespXY <- read.csv2("myRespXY_corrected.csv",header=TRUE, sep=",",stringsAsFactors = F)

# TO CONSIDER: if you want to set some data aside
# Use eval.resp.var, eval.resp.xy and eval.expl.var

# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, # Presence-Absence vector
                                     expl.var = myExpl,
                                     resp.xy = myRespXY, # Data frame longitude latitude (must be same amount of rows as resp.var)
                                     resp.name = myRespName,
                                     filter.raster = TRUE) #Try also FALSE if there is no over representation of one location


myBiomodData
plot(myBiomodData)
table(myBiomodData@data.species)
show(myBiomodData)
str(myResp)



# Transform true absences into potential pseudo-absences
myResp.PA <- ifelse(myResp == 1, 1, NA)  # Remove true absences


# Format Data with pseudo-absences : random method
myBiomodData.PA <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                      expl.var = myExpl,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName,
                                      PA.nb.rep = 4, # Generate 4 pseudo-absence data sets
                                      PA.nb.absences = 10000, # Number of pseudo-absences
                                      PA.dist.min = NULL,  # Minimum distance between PA and presence
                                      PA.strategy = 'random', # "random" or "disk" (avoid placing pseudo-absences near presences)
                                      filter.raster = TRUE)

myBiomodData.PA
plot(myBiomodData.PA)
# Check pseudo-absence count
table(myBiomodData.PA@data.species)












# --------- Cross validation ------
# # k-fold selection
# cv.k <- bm_CrossValidation(bm.format = myBiomodData,
#                            strategy = "kfold",
#                            nb.rep = 2,
#                            k = 3)

# cv.r <- bm_CrossValidation(bm.format = myBiomodData,
#                            strategy = "random",
#                            nb.rep = 5,
#                            perc = 0.7) # keep 70% of the data for evalu

# 
# # stratified selection (geographic)
# cv.s <- bm_CrossValidation(bm.format = myBiomodData,
#                            strategy = "strat",
#                            k = 2,
#                            balance = "presences",
#                            strat = "x")
# head(cv.k)
# head(cv.s)

# ------ Retrieve modeling options ----

# # big boss parameters
# opt.b <- bm_ModelingOptions(data.type = 'binary',
#                             models = c('SRE', 'XGBOOST'),
#                             strategy = 'bigboss')
# 
# # tuned parameters with formated data
# opt.t <- bm_ModelingOptions(data.type = 'binary',
#                             models = c('SRE', 'XGBOOST'),
#                             strategy = 'tuned',
#                             bm.format = myBiomodData)
# 
# opt.b
# opt.t

# changing some options for the models to stop over fitting
# then add bm.options = myBiomodOption in BIOMOD_Modeling()
# myBiomodOption <- BIOMOD_ModelingOptions(GAM = list( algo = 'GAM_mgcv',
                                                  #   type = 's_smoother',
                                                   #  k = 4),
                                        # MAXENT = list(
                                         #  maximumiterations = 500,
                                          # threshold = FALSE))

# --------- Run modelisation -------

# Model single models
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format    = myBiomodData,
  models       = c('GLM', 'GAM', 'RF', 'GBM'),  # Full list: ANN, CTA, FDA, GAM, GBM, GLM, MARS, MAXENT, MAXNET, RF, SRE, XGBOOST
  modeling.id  = 'GLM&GAM&RF&GBM_world',
  CV.strategy  = 'random', # type of CV, Random partitioning into training and test sets (recommended). can also be kfold (small datasets), LOO (slow, small ds), block (spatial bias and autocorrelation)
  CV.nb.rep    = 5, # btw 3 and 5 good balance for medium datasets, 10: more stability but more computation
  CV.perc      = 0.8, # btw 0.7 - 0.8: common practice (70-80% training, 20-30% testing)
  OPT.strategy = 'bigboss', # bigboss optimizes hyper parameters, "grid_search" performs a grid search for best parameters (accurate but slow), user.defined is manually set params 
  var.import   = 5, # Number of permutations to assess variables importance (higher is more precise but slower. 3-5 for best balance)
  metric.eval  = c('TSS','ROC'))

# seed.val = 123)
# nb.cpu = 8)
myBiomodModelOut


# Extract individual models
# Load individual models from the BIOMOD2 output
#rf_models <- BIOMOD_LoadModels(myBiomodModelOut, algo = "RF")

# Get evaluation scores & variables importance
eval <- get_evaluations(myBiomodModelOut) # 
varimp <- get_variables_importance(myBiomodModelOut)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)], # Here you can change the subset of models chosen in c()
                      fixed.var = 'median') # When generating the response curve for a given variable, all other variables are set to their median value.
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[3],
                      fixed.var = 'median',
                      do.bivariate = TRUE)



# ---- if satisfied with models: ENSEMBLE!! ----

myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean'), #'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.5),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3)
                                      #EMci.alpha = 0.05,
                                      #EMwmean.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = "all",
                      fixed.var = 'median',
                      do.bivariate = TRUE)
