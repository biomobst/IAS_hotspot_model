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
library(glue)
library(readxl)# to rename files
setwd("/home/gunnarandersson/Dokument/Projekt/HAV2025") 
# Explore datasets in the package
#list_datasets()

# Explore layers in a dataset
#list_layers() 


projectpath <- "~/Dokument/Projekt/HAV2025"
dir <- "~/Dokument/Projekt/HAV2025/data/Biooracle.download/"

#
### BIO-ORACLE LAYERS
# Explore environmental variables from Bio-Oracle
layers.bio2 <- as.data.frame(list_layers(simplify=F) )#simplify=F if need more info
layers.bio2$dataset_id

#remove irrelevant
exclude <- c(grep("terrain",layers.bio2$dataset_id),grep("kdpar",layers.bio2$dataset_id))
layers.bio2<- layers.bio2[-exclude,]

#Find all variables
dataset_variables <- unique(lapply(layers.bio2$dataset_id, function(i)
  strsplit(i,"_")[[1]][1]
))
dataset_scenarios <- unique(lapply(layers.bio2$dataset_id, function(i)
  strsplit(i,"_")[[1]][2]
))
dataset_scenarios<- dataset_scenarios[-8]#remove mean
dataset_scenarios_titles <- unique(layers.bio2$title)

#get variables from par
info_layer("par_mean_baseline_2000_2020_depthsurf")



#layers.bio2$dataset_id[grep("_mean",layers.bio2$dataset_id)]

#get variables from layer
#layers.bio2[1,]
#layers.bio2[1,c("title")]

#testlayer <- layers.bio2$dataset_id[1]
#testlayerspellout <- layers.bio2$title[1]


###
#scenario <- dataset_scenarios[1]
#selected.datasets.nr <-grep(scenario,layers.bio2$dataset_id )
#layers.bio2$dataset_id[selected.datasets.nr ]
#for(i in selected.datasets.nr){
#  print(paste(layers.bio2$title[i], layers.bio2$dataset_id[i], sep=";"))
 # print(paste(unlist(info_layer(layers.bio2$dataset_id[i])[["variables"]]["variable_name"]),collapse=";"))
#}


variable.selection <- read_excel("data/dataset och variabler 2025.xlsx")
#remove not used
variable.selection <- variable.selection[!is.na(variable.selection$Final),]
variable.selection <-variable.selection[-3,]# chl does not exist at depthmean ssp119

for(sel.sen in c(4:7)){
scenario <- dataset_scenarios[[sel.sen]][1]

ssp_scenarios <- list()
constraints.base <- list(
  #latitude = c(25, 80),
  #longitude = c(-15, 40),
  time = c("2010-01-01T00:00:00Z", "2010-01-01T00:00:00Z") # Ensure time is correctly specified
)
constraints.ssp <- list(
  #latitude = c(25, 80),
  #longitude = c(-15, 40),
  time = c("2020-01-01T00:00:00Z", "2020-01-01T00:00:00Z") # Ensure time is correctly specified
)


i <-1 #
#variable.selection$Spellout

if(scenario == "baseline"){constraints <- constraints.base}else{constraints <- constraints.ssp}
scenario.info <- list()
scenario.info[[scenario]]<-list()
scenario.info[[scenario]][["datasets"]]<- list()
#remove par for otehr than baseline
number = ifelse(scenario == "baseline",length(variable.selection$Spellout),length(variable.selection$Spellout)-1)
for(i in 1:number){
scenario.info[[scenario]][["datasets"]][[i]]<- list(dataset_id = variable.selection$dataset[i] ,
                                                    variables = gsub(";",",",variable.selection$Final[i]),
                                                    constraints = constraints)
}

for (dataset.nr in 1:length(scenario.info[[scenario]][["datasets"]])) {
#for (dataset.nr in 1:3) {
  
  dataset <- scenario.info[[scenario]][["datasets"]][[dataset.nr]]
  dataset_id <- dataset$dataset_id
  if(!(scenario == "baseline")){
    scenario.pos <- which(!is.na(match( strsplit(dataset_id,"_")[[1]],"baseline")))
  variable.vector <- strsplit(dataset_id,"_")[[1]]
  variable.vector[scenario.pos]<- scenario
  variable.vector[scenario.pos+1]<- 2020
  variable.vector[scenario.pos+2]<- 2100
  
  dataset_id <- paste(variable.vector, collapse ="_")
  }#end if
  strsplit(dataset_id,"_")[[1]]
  variables <- dataset$variables
  constraints <- dataset$constraints

  for(variable in strsplit(variables,",")[[1]] ){
  # Download and load data layers of choice 
 # dir <- "/home/gunnarandersson/Dokument/Projekt/HAV2025/data/Biooracle.download/"
  if (!dir.exists(dir)) dir.create(dir) 
  
  info_layer(dataset_id)
  a <- download_layers(dataset_id, variables = variable, constraints = constraints, directory= dir)
  filename_with_ext <- basename(terra::sources(a))[1]
  file.rename(
    from = glue("{dir}{filename_with_ext}"),
    to = glue("{dir}{dataset_id}_{variables}.nc"))
  b <- brick(glue("{dir}{dataset_id}_{variables}.nc"))
  filename1 <- glue("{dir}{dataset_id}_{variables}.nc")
  
  # datalayer.tiff for tif, etc
  outdir <- paste(dir,"datalayer.nc/",scenario,"/",sep="")
  if (!dir.exists(outdir)) dir.create(outdir) 
  depth <- strsplit(dataset_id,"_")[[1]][length(strsplit(dataset_id,"_")[[1]])]
  filename <-paste(outdir,variable,"_",depth,".tif",sep="")
  # to save as tiff
  #filename <-paste(outdir,variable,"_",depth,".tif",sep="")
  #writeRaster(b, filename, bylayer=TRUE)
  #to keep as nc
  filename <-paste(outdir,variable,"_",depth,".nc",sep="")
  file.rename(filename1,filename)
  }#for var
}
}#end sel.sen

## rename files
for(sel.sen in c(2:7)){
  #sel.sen <- 1
  scenario <- dataset_scenarios[[sel.sen]][1]
  outdir <- paste(dir,"/datalayer.tiff/",scenario,"/",sep="")
  a = list.files(outdir, pattern = "*_1.tif")
  for(i in a){
    i <- paste(outdir,i,sep="")
    file.rename(i,gsub( "_1.tif", ".tif",i))
  }
}