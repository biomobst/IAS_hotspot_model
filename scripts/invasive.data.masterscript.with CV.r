#invasive.data.masterscript.r
source("SEanalytics.functions.r")

require(raster)
require(rgdal)
require(zip)

########################################################################################
### Ecological niche modelling for invasive species
### Gunnar Andersson, National Veterinary Institute, SVA, Sweden  gunnar.andersson@sva.se
### Prepared fro SE-analytics Mars 2020
########################################################################################


stringsAsFactors = F
#Folder with the original rasterdata
Rasterpaths <- c("C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Benthic_Avg/Present_Benthic_Avg",
                 "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Surface/Present_Surface")

#Folder where the rasterstacks are stored
Stackpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/rasterstacks"

# Folder containg presence and absence data for each species and preudoabsence data used for all species.
#Speciespath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/species to use apr 2020"
Speciespath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/species to use"

#This folder contains files describing the cross validation scheme for each species
#Iterations.path <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/iterations/iterations apr 2020"
Iterations.path <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/iterations"

ROC.path <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ROC"
#ROC.path <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ROC/roc apr 2020"



# where plots should be stored
Plotpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/plots" 
#Plotpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/plots/plots apr 2020" 


#The "outpath" contains the .csv files with species data and extracted environmental variables
Outpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/RF.indata/"

# the modelpath is where the models from random forest is saved. The files will also contain predictions from the cross validation,
# used as input when ROC curves are prpares
Modelpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/models"

#Modelpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/models/models apr 2020"

resultpath <-  "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/Resultat"

# This contains the predicted species distribuition maps as rda files. These are used input for the plots.
mappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/maps"
#mappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/maps/maps apr 2020"

# stores the prediction resutls as geotiff rasters. Both data for individual species and for the different measures
rastermappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters"
#rastermappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters/predict rasters apr 2020"

#Paths to folders where shapefiels are stored describing ICES areas. Used as overlays for plots only.
ICESpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ICES_areas"
ICESecopath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ICES_ecoregions"
ICESstatrectpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ICES_StatRec_mapto_ICES_Areas"

# read in the region map used in the plots, and transfor to WGS84 - World Geodetic System 1984
# actually the map will not be needed until later. Readin git could be omitted here to save memory.
shape2 <- readOGR(dsn="H:/R gruppen/Kartor/NUTS_shapefile",layer="nutsByHand")
proj4string(shape2)#GRS80 EPSG:4019
shape2 <- spTransform(shape2, CRS("+init=epsg:4326"))
plot(shape2)# just to check that it appears right
#########################################################

# Read the
#Data.table <- read.csv2("data.table.csv",header=TRUE, sep=";",stringsAsFactors = F)
#Data.table <- read.csv2("data.table.apr2020.csv",header=TRUE, sep=";",stringsAsFactors = F)
Data.table <- read.csv2("data.table.Anteaeolidiella.lurana2020.csv",header=TRUE, sep=";",stringsAsFactors = F)


names(Data.table)[1]<- "species" # just to check

# select one species to try out the code. Not used in the loop.
Species = Data.table$species[1] #"Ficopomatus enigmaticus"
#Define which stack to used when extracting environmental data-
Stack <- "globalStack.rda"#"globalStack.rda" or "europeStack.rda"

#Prevent that strings in csv files are red in as factors variables
stringsAsFactors= FALSE

#for(Species in Data.table$species[1:2]){
  for(Species in Data.table$species){
    
  print(Species)
  # Read in present absent and pseudoabsent points, 
    #convert these points to spatial coordinates and extract environmental variables from rasterstack
 species.data <-  read.and.extract(data.table = Data.table,
                                   species = Species,
                                   stack = Stack, 
                                   speciespath = Speciespath,
                                   stackpath = Stackpath,
                                   plotpath = Plotpath,
                                   outpath = Outpath)
 
 ### Check that there are no erroneous entries in the files. 
 ## the most common synonyms for "present" and "absent" are identified and entries harmonized.
 #For entires where status cannot be determiend the line is removed
 # If data was checked in the earleir stage these lines should not find any mistakes
 my.data <- species.data
 print(paste(Species,paste(unique(my.data$occurrenceStatus) ))  )
 my.data$occurrenceStatus <- as.character(my.data$occurrenceStatus)
 present.synonyms <-  which(my.data$occurrenceStatus == "present"|my.data$occurrenceStatus == "Present"| my.data$occurrenceStatus == "established"| my.data$occurrenceStatus == "Established")
 my.data$occurrenceStatus[present.synonyms] <- "present"
 absent.synonyms <-  which(my.data$occurrenceStatus == "Absent"| my.data$occurrenceStatus == "")
 
 my.data$occurrenceStatus[absent.synonyms] <- "absent"
 remove <- which(!my.data$occurrenceStatus == "absent"  & !my.data$occurrenceStatus == "present")
 print(paste("remove",remove))# promt if lines are removed
 if(length(remove > 0)){
   my.data <- my.data[-remove,]                
 }
 species.data <- my.data
 # Save the species data with environmental variables in selected folder
  write.csv2(species.data, file = paste(Outpath,"/",Species,"indata.csv", sep=""),row.names=F)
}
  #print(species)
###

##################################################################
### prepare iterations
##################################################################


#Call the split data function to generate a list object that describes which entries are 
# used as training data in each repetition and CV-fold in the cross-validation
# the resutls is stored as a rda file in the folder defined in "iterations path"

#The settings are made in teh "split data" function
#  n.repeat <- 5 this defines the number of times the cross validation process is performed
#  CV.level <- 5 this defienes the level of the cross validation
# to speed up thing it is possible to make  2-fold cross validation and just one repeat.
# when the data is split into ICES statistical rectangles based on their latitude and longitide.
#If data from other parts of the globe are used this may have to be modified
#

#for(Species in Data.table$species[1]){
for(Species in Data.table$species){
    
  # Species <- Data.table$species[1] 
 # iterations.path
  split.data(species = Species, 
             indata.path = Outpath,
             iterations.path= Iterations.path)

}
####################
##################################################################
  ## Run random forests. the function will prepare data and execute.
##################################################################
for(Species in Data.table$species){
  #Species <-  Data.table$species[3]
#for(Species in Data.table$species[1:2]){
  # the object "rf.output.list contains all resutls from random forests as a list object
  rf.output.list <- run.random.forests(species = Species, 
                selvar = "all",
                indata.path = Outpath,
                iterations.path= Iterations.path)
  
  # run a garbage collection to free memory
  gc()
  # the rf.output.cv object contains the results from the cross validations analysis. 
  #Predictions for each entrie are stored for each repeat and cv-fold. The random forest model is not saved, to save memory
  rf.output.cv <- rf.output.list[["RF.results.CV"]]
  # this object contains the random forests model obtainied when all indata is used as trainingset.
  # the predicted probability of being present is stored for each data poing
  rf.output <- rf.output.list[["RF.results.alldata"]]
  
  #print results to see progression
  names(rf.output)
  print(rf.output$"RF.selected")
 # print(rf.output$"response.sel")
  
  #save cross validation results for later calculation of ROC
  save(rf.output.cv, file= paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",Species,".rda",sep=""))
  #save the random forest model for late prediction of maps
  save(rf.output, file= paste(Modelpath,"/RF.model.and.predictions.eur.wt.",Species,".rda",sep=""))
  #remove the large object generated to make space for the next species
  rm(rf.output.list)
  rm(rf.output.cv)
  rm(rf.output)
  gc()
  
}
##################################################################
# calculate ROC curves and plot 
##################################################################
for(Species in Data.table$species){
  #Species <-  Data.table$species[2]
  load(paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",Species,".rda",sep=""))
  load(paste(Iterations.path,"/occurance.iters_",Species,".rda",sep=""))
  
  lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
  my.data <- read.csv2(lista.csv[grep(Species,lista.csv)],header=T)
  
  .true.class <-my.data$occurrenceStatus
  
  #######recreate the experiment setup
  nrep <- length(all.occurance.iters)
  CV.level <- length(all.occurance.iters[[1]])
  process <- seq(1,nrep*CV.level)
  rep <- sort(rep(seq(1,nrep),CV.level))
  iter <- rep(seq(1,CV.level),nrep)
  process.plan <- cbind(process,rep,iter)# may be obsolete when not running parallell script
  ###
  all.rep <- unique(process.plan[,2])
  curr.rep <- 2
  i <-1
   method <- "allvars"
  #### calculate one ROC curve per repetitio of the Cross validation
   all.ROC <- lapply(all.rep, function(curr.rep)
     calc.ROC(   rf.output.cv[[curr.rep]],.true.class)
   )
  
   # prepared a average ROC curve based on all repetitions togeather
    mean.ROC <- calc.ROC(lapply(1:length(process.plan[,1]), function(i)
      rf.output.cv[[   process.plan[i,2]  ]][[ process.plan[i,3]  ]] ) ,.true.class)
   
    # .RF.result <-lapply(1:length(process.plan[,1]), function(i)  rf.output.cv[[   process.plan[i,2]  ]][[ process.plan[i,3]  ]] )
   #      dim( rf.output.cv[[curr.rep]][[i]]$predict.selected )
#rf.output.cv[[curr.rep]][[i]]$predict.selected[1,]
  
    #Plot the ROC cuve in the dedicated folder
  plot.ROC(.ROC.path=ROC.path,
           .my.species = Species,
           .mean.ROC = mean.ROC,
           .all.ROC =all.ROC)
  #remove objects t make list for next run
    rm(all.occurance.iters)
    rm(rf.output.cv)
    rm(all.ROC)
    rm(mean.ROC)
    gc()
}
##################################################################
## plot maps
##################################################################
#plots
#load rasterlayers
load(paste(Stackpath,"europeStack.rda", sep="/") )
plot(var2[[1]])
bar <- stack(var2)
colorsBrBG <- rev(divPalette(n=12, name = c( "BrBG") ) )
colorsBrBG2 <- rev(divPalette(n=18, name = c( "BrBG") ) )

plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[10:18])
colorsBlue <- seqPalette(n=,12, name = c( "Blues") ) 
shape.ecoregions <- readOGR(dsn=ICESecopath,
                            layer="ICES_ecoregions_20171207_erase_ESRI",encoding="UTF-8")
brk <- c(seq(0, 1,by=0.1),1.05)


# Species  <- Data.table$species[1]

#call the function to predict the species distribution at raster level.
#The function return the prediction as a raster object but also make plots as .png 
# the lines between png() and dev.off() may be removed/inactivated if the png plots are not wanted.
for(Species in Data.table$species){
 map <-  plot.maps(species = Species,
            indata.path = Outpath,
            modelpath = Modelpath,
            plotpath = Plotpath, 
            colors = c(colorsBr, "lightgrey"),
            brk = brk
            )
 save(map, file= paste(mappath,"/predicted.map.",Species,".rda",sep=""))
}

#plot(map)
### compare maps
#Species  <- Data.table$species[1]
colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )
colorsBrBG2 <- rev(divPalette(n=22, name = c( "BrBG") ) )

#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[12:22], "lightgrey")


#brk <- c(seq(0, 1,by=0.1),1.1)

brk.log <- c(seq(-3, 0,by=0.25),0.25)
length(brk.log)

# define the area to plot. In this case the coordinates fro teh swedish map
xlim=c(0,30)
ylim = c(50,70)
rm(map)
load(paste(Stackpath,"europeStack.rda", sep="/") )
plot(var2[[1]])
bar <- stack(var2)

for(Species in Data.table$species){
  #reload rf.output
   load( paste(mappath,"/predicted.map.",Species,".rda",sep="") )
  #load(paste(modelpath,"/RF.model.and.predictions.eur.wt.",species,".rda",sep=""))
 
 #map[is.na(mean(bar))]<- 0.000001
  lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
  my.data <- read.csv2(lista.csv[grep(Species,lista.csv)],header=T)
  
  # prepare a rasterlayer with the 10-log of the predicted probability of presence.
  # An arbitrary small number is added, as log10(0) is not defined
  logmap <- log10(map+0.001) #log10(max(map , 0.001, na.rm=T))
  
  # an arbitrary large number is inserted at locations where presence is undefined, in other words land.
 logmap[is.na(mean(bar))]<- 0.24
 # save(logmap, file= paste(mappath,"/predicted.map.log.",Species,".rda",sep="") )
  
 # map[is.na(map)]<- -0.1
 
 # plot maps at Eurpean and Swedish scale this presence points added.
 # Plot without indicating xlim and ylim. This gives a plot area defined by the raster extent
  png(paste(Plotpath,"/",Species, ".logprob.Eur.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(logmap,col=colorsBr, breaks=brk.log)
  plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=0.3, lwd=0.4)
  dev.off()
  ### same for Sweden
  # make the plot using the xlim and ylim defined earlier.
  png(paste(Plotpath,"/",Species, ".logprob.Swe.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(logmap,col=colorsBr, breaks=brk.log , xlim=xlim, ylim=ylim )
  plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=0.3, lwd=0.25)
  dev.off()
  
  
  
  # save predicted probabilities as rasterfiles.
  newfile <- paste(rastermappath,"/linear.prob.",Species,sep="")
  #print(paste(file, newfile , sep="###########"))
  writeRaster(map, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
  
  newfile <- paste(rastermappath,"/log.prob.",Species,sep="")
  #print(paste(file, newfile , sep="###########"))
  writeRaster(logmap, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
  
  rm(logmap)
  rm(map)
  gc()
}

############# combine layers
# define color palettes. White and lightgrey to incidate values outside range. 
colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )
colorsBrBG2 <- rev(divPalette(n=22, name = c( "BrBG") ) )

#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[12:22], "lightgrey")


#brk <- c(seq(0, 1,by=0.1),1.1)

# the cutoffs to used in the color scale at log-scale
brk.log <- c(seq(-3, 0,by=0.25),0.25)

#Find the predicted maps for the different species
lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.*",sep="/"))

#  read in some extra maps, from an older run
#  rastermappath2<- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters"
# lista.ras2 <- Sys.glob(paste(rastermappath2,"*linear.prob.*",sep="/"))
# lista.ras <- c(lista.ras, lista.ras2)
#make a rasterstack of all predicted species, fur different computations

predvar<-stack(lista.ras)
# calculate the mean probability of presence for all species, as a rasterlayer.
mean.predvar <- mean(predvar, na.rm=T)

newfile <- paste(rastermappath,"/mean.prob.all.species",sep="")
writeRaster(mean.predvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)

logpredvar <- log10(mean(predvar, na.rm=T) + 0.001)

newfile <- paste(rastermappath,"/log 10 of mean.prob.all.species",sep="")
writeRaster(logpredvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)

logpredvar[is.na(mean(bar))]<- 0.25
png(paste(Plotpath,"/mean.all.eur", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(logpredvar,col=colorsBr, breaks=brk.log)
plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()
png(paste(Plotpath,"/mean.all.Swe", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(logpredvar,col=colorsBr, breaks=brk.log,xlim=xlim, ylim=ylim)
plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()

### mean of log10 predvar
custom.pred <- log10(predvar[[1]]+0.001)
for(i in 2:length(names(predvar))){
  custom.pred <- custom.pred+ log10(predvar[[i]]+0.001)
}
plot(custom.pred)
newfile <- paste(rastermappath,"/sum of log10 prob all.species",sep="")
writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
custom.pred <- custom.pred/length(names(predvar))
plot(custom.pred)
newfile <- paste(rastermappath,"/mean of log10 prob all.species",sep="")
writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)


#######################################
# customzed combined layers
#Dessa kartor är nog viktigast. Jag undrar om vi kan få rasterfiler för
#2. medelvärdet av Log-sannolikheten
#3. summerad Log-sannolikheten
#4. antal arter över 50%-tröskel delad upp efter 2 grupperingar (group:Environment, group:Salinity) enligt bifogat masterfile (IASmodellng_speciesList)
#
combine.table <-  read.csv2("IASmodelling_speciesList_20200107.csv",header=TRUE, sep=";",stringsAsFactors = F)
names(combine.table)<- c("Vetenskapligt namn",	"Svenskt namn",	"Source",
                          "Samlat riskutfall"	,"included"	,"Environment","Taxonomy",	"Salinity"	,"Habitat",	"gbif occurence",	"useful points",	"comments")

lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.*",sep="/"))



lista.ras.copy <- gsub("_"," ", lista.ras)# to be able to match species names and sort list

in.ras.order <- sapply(1:length(Data.table$species), function(i) which(grepl(Data.table$species[i] ,lista.ras.copy)))


# layers soted in the order they appear in the Data table"
lista.ras.sort <- lista.ras[in.ras.order]
lista.ras.copy.sort <- lista.ras.copy[in.ras.order]


# sort the information on how to combine in the same order as the other tables
in.modelling.order <- sapply(1:length(Data.table$species), function(i) which(grepl(Data.table$species[i] ,combine.table$"Vetenskapligt namn")))
combine.table.sort <- combine.table[in.modelling.order,]


# extract species names from lista ras
# lista.ras.sort <- lista.ras
# lista.ras.sort.names <- sapply(1:length(lista.ras), function(i) strsplit(lista.ras[i], "prob[.]|[.]tif")[[1]][2])
#in.modelling.order <- sapply(1:length(lista.ras.sort.names), function(i) which(grepl(lista.ras.sort.names[i] ,combine.table$"Vetenskapligt namn")))
# combine.table.sort <- combine.table[in.modelling.order,]


# prepare sorted stack
#unique(combine.table.sort$"Salinity")
#define what species are included in what stack
sel <- list()
sel[["all species"]] <- rep(1,length(combine.table.sort$"Environment"))
sel[["freshwater"]] <- ifelse(combine.table.sort$"Salinity" == "estuarine/freshwater" ,1,0)
sel[["marine environment"]] <- ifelse(combine.table.sort$"Salinity" == "marine" ,1,0)
for(group in unique(combine.table.sort$"Environment") ){
  sel[[as.character(group)]] <- ifelse(combine.table.sort$"Environment" == group ,1,0)
}
sel[["zoo planct + bent"]] <- sel[["phytobenthos"]] + sel[["zooplankton"]]
sel[["phyto planct + bent"]] <- sel[["phytoplankton"]] + sel[["phytobenthos"]]
sel[["plankton zoo + phyt"]] <- sel[["phytoplankton"]] + sel[["zooplankton"]]
sel[["benthos zoo + phyt"]] <- sel[["phytobenthos"]] + sel[["zoobenthos"]]


sortpredvar <-stack(lista.ras.sort)
#plot(mean(sortpredvar))
### weighet average
#logpredvar <- log10(mean(predvar, na.rm=T) + 0.001)


combine.table.sort$"Samlat riskutfall"[which(is.na(combine.table.sort$"Samlat riskutfall"))] <- 0
combine.table.sort$"Samlat riskutfall"[is.na(combine.table.sort$"Samlat riskutfall")] <- "0"
combine.table.sort$"Samlat riskutfall"[which(combine.table.sort$"Samlat riskutfall"== "unknown")] <- "2.5"
combine.table.sort$"Samlat riskutfall"[which(combine.table.sort$"Samlat riskutfall"== "")] <- "0"

weighting <- c("weighted","raw")

for(mode in weighting){
for( g in 1:length(names(sel)) ){
  group <- names(sel)[g]
  if(mode == "weighted"){
  weight.map <-(3 + log10(sortpredvar[[1]]+ 0.001) ) * as.numeric(combine.table.sort$"Samlat riskutfall"[1])*sel[[group]][1]
  }else{
    weight.map <-(3 + log10(sortpredvar[[1]]+ 0.001) ) * 1*sel[[group]][1] 
  }
  
  for(i in 2:length(in.modelling.order)){
    if(mode == "weighted"){
  weight.map <- weight.map +( (3+log10(sortpredvar[[i]]+ 0.001)) * as.numeric(combine.table.sort$"Samlat riskutfall"[i]) *sel[[group]][i]) 
    }else{
      weight.map <- weight.map +( (3+log10(sortpredvar[[i]]+ 0.001)) *sel[[group]][i]) 
      
    }
  }
  if(mode == "weighted"){
weight.map <- weight.map/sum(as.numeric(combine.table.sort$"Samlat riskutfall")[which(sel[[group]]>0.5)])
  }else{
    weight.map <- weight.map/length(which(sel[[group]]>0.5) )
    
}
#print(paste(group,"weights",as.numeric(combine.table.sort$"Samlat riskutfall")[which(sel[[group]]>0.5)]))
#weight.map[is.na(weight.map)]<- -1
#plot(weight.map)

newfile <- paste(rastermappath,"/",mode,".log.map.",group,sep="")
print(newfile)
#print(paste(file, newfile , sep="###########"))
writeRaster(weight.map, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
####
}#end all groups
} #end mode
#1. antal arter över 50%-tröskel


####################################################################
# maps with % treshold
for(cutoff in c(0.10,0.33,0.50)){
for(mode in weighting){
  for( g in 1:length(names(sel)) ){
    group <- names(sel)[g]
    
  #sel <- ifelse(combine.table.sort$"Environment" == group ,1,0)
  
    if(mode == "weighted"){
      sum.50.map.sel <- ceiling(sortpredvar[[1]]-cutoff )*sel[[group]][1] *
        as.numeric(combine.table.sort$"Samlat riskutfall"[i])
      }else{
    sum.50.map.sel <- ceiling(sortpredvar[[1]]-cutoff )*sel[[group]][1]  #*as.numeric(combine.table.sort$"Samlat riskutfall"[1])
      }
     for(i in 2:length(in.modelling.order)){
       if(mode == "weighted"){
         sum.50.map.sel <- sum.50.map.sel+(ceiling(sortpredvar[[i]]-cutoff )*sel[[group]][i] )* 
           as.numeric(combine.table.sort$"Samlat riskutfall"[i])
       }else{
    sum.50.map.sel <- sum.50.map.sel+(ceiling(sortpredvar[[i]]-cutoff )* sel[[group]][i] )# *as.numeric(combine.table.sort$"Samlat riskutfall"[i])) 
       }
    }
  #plot(sum.50.map.sel)
  newfile <- paste(rastermappath,"/sum species ",mode," ",cutoff," ",group,sep="")
  writeRaster(sum.50.map.sel, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
  
  }
}
}
################################## plot selected layers

colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )



groups <- c("*sum species*","*weighted.log*","*raw.log*")

for(group in groups){
#brk <- c(seq(0, 1,by=0.1),1.1)
lista.ras <- Sys.glob(paste(rastermappath,group,sep="/"))
plotvar <- stack(lista.ras)

for(lager in 1:length(lista.ras)){

#lager <- 1
#plot(plotvar[[lager]])

min <- min(values(plotvar[[lager]]), na.rm =T)
max <- max(values(plotvar[[lager]]), na.rm =T)

nval <- ceiling( max - min )

plotvar[[lager]][is.na(plotvar[[lager]])] <- max+(max-min)/40
colorsBrBG2 <- rev(divPalette(n=2*nval, name = c( "BrBG") ) )

#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[seq((nval+1),2*nval)], "lightgrey")

#length(colorsBrBG2[12:22])

brk.custom <- c(seq(min, max,length.out= length(colorsBr)) ,max+(max-min)/20)

content1 <- strsplit(lista.ras[lager],"/")
content2 <- content1[[1]][length(content1[[1]])]
content3 <- strsplit(content2,"[.]tif")[[1]][1]

newfile <- paste(Plotpath,"/XX",content3,".EUR.png",sep="")
png(newfile,  width = 180, height = 180, units = "mm", res=1200)

# define legend
arg <- list(at=c(min, max), labels=c(min,max))
plot(plotvar[[lager]],col=colorsBr, breaks=brk.custom,axis.args=arg)
plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
dev.off()

newfile <- paste(Plotpath,"/XX",content3,".SWE.png",sep="")
png(newfile,  width = 180, height = 180, units = "mm", res=1200)
plot(plotvar[[lager]],col=colorsBr, breaks=brk.custom,axis.args=arg,xlim=xlim, ylim=ylim)
plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()
}#end lager
}#end gropus

################################################################################
#                              END                                         #####
#¤#################################################################################
### pca on rasters
#https://www.rdocumentation.org/packages/RStoolbox/versions/0.2.6/topics/rasterPCA
library(ggplot2)
library(reshape2)
library(RStoolbox)
#Just a test to run PCA on the rastetrlayers... 
## Run PCA
#rasterPCA(img, nSamples = NULL, nComp = nlayers(img), spca = FALSE, maskCheck = TRUE, ...)
set.seed(25)
rpc <- rasterPCA(predvar,nComp = 3,spca= T)
summary(rpc$model)
loadings(rpc$model)


ggRGB(rpc$map,1,2,3, stretch="lin", q=0)
  plots <- lapply(1:3, function(x) ggR(rpc$map, x, geom_raster = TRUE))
 plot(plots[[1]],xlim=xlim, ylim=ylim)
 plot(plots[[2]],xlim=xlim, ylim=ylim) 
 plot(plots[[3]],xlim=xlim, ylim=ylim)
####### obtain oub prediction performance
 for(Species in Data.table$species){
   file <-  paste(Modelpath,"/RF.model.and.predictions.eur.wt.",Species,".rda",sep="")
   #print(file)
   load(file)
   print("########################")
   print("")
   print(Species)
   print(rf.output$"RF.selected"$"confusion")
   
 }
