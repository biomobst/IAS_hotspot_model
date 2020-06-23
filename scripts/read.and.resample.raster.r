#extract.data.r
require(raster)
require(rgdal)
#install.packages("zip")
require(zip)

#This script reads zipped rasterlayers and prepares a rasterstack, for use in subsequent analysis

#Paths to the different sets of layers
rasterpaths <- c("C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Benthic_Avg/Present_Benthic_Avg",
                 "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Surface/Present_Surface")


# Save the rasterstack in this directory
stackpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/rasterstacks"

#############################
# read raster data
dir.create(tmp <- tempfile())
dir.create(tmp2 <- tempfile())

#find all layers in the different folders. Copy to temporary directory for processing.
for(mypath in rasterpaths){
  lista.ras<- Sys.glob(paste(mypath,"*.tif.zip",sep="/"))
  for(file in lista.ras){
    unzip(file,exdir =tmp)
  }
  print(lista.ras)
}



lista.ras2 <- Sys.glob(paste(tmp,"*.tif",sep="/"))

#Compare the rasterlayers to see that they have the same projections, and thus can be stacked. Exclude layers that do not comply
outcome <- c()
for(i in 1:length(lista.ras2)){
  outcome <- c(outcome, try(compareRaster(raster(lista.ras2[1]), raster(lista.ras2[i]))) )
}

# Prepare a stack of the layers that are in same format as the first one
var<-stack(lista.ras2[ which(outcome == "TRUE")])
plot(mean(var))


# read in all data so that it will be saved in the .rda file
var <- readAll(var)

# Save the stack of functioning rasters on global map
save(var, file=paste(stackpath,"globalStack.rda", sep="/"))

#list any non functional layers
lista.ras2[  which(!outcome == "TRUE") ]

# Just a check of non compilant rasters
compareRaster(raster(lista.ras2[  which(!outcome == "TRUE") ][1]), raster(lista.ras2[  which(!outcome == "TRUE") ][2]) )

#get extent of rasters 
extent(raster(lista.ras2[  which(!outcome == "TRUE") ][1]))# raster with non standard projection
extent(raster(lista.ras2[  which(outcome == "TRUE") ][1]))#rasters with standard prjection
#plot(raster(lista.ras2[  which(outcome == "TRUE") ][1]))

#get projectins for rasters
proj4string(raster(lista.ras2[  which(!outcome == "TRUE") ][1]))#non - standard, not used
proj4string(raster(lista.ras2[  which(outcome == "TRUE") ][1]))#standard
############# beskär och omsampla  

#Define extent of the map to be used. In 
e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax 
lista.ras.tmp <- Sys.glob(paste(tmp,"*.tif",sep="/"))
for(file in lista.ras.tmp[-c(91:92)]){
  copen<-raster(file)
  #proj4string(copen )
  #extent(copen)
  rc <- crop(copen, e)	
  tmp.string <- strsplit(as.character(tmp), "[\\]")[[1]][length(strsplit(as.character(tmp), "[\\]")[[1]])]
  tmp.string2 <- strsplit(as.character(tmp2), "[\\]")[[1]][length(strsplit(as.character(tmp), "[\\]")[[1]])]
  
  newfile <- gsub( tmp.string,tmp.string2, file)
  #print(paste(file, newfile , sep="###########"))
  writeRaster(rc, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
  
}

lista.ras.tmp2 <- Sys.glob(paste(tmp2,"*.tif",sep="/")) #Sys.glob(paste(tmp2,"*.*",sep="/"))





####
require(fBasics)
colorsBlue <- seqPalette(n=,10, name = c( "Blues") ) 
colorsBrBG <- rev(divPalette(n=10, name = c( "BrBG") ) )

#plot(raster(lista.ras.tmp2[  which(!outcome == "TRUE") ][1]), col=colorsBlue)
plot(raster(lista.ras.tmp2[1]), col=colorsBrBG)

var2<-stack(lista.ras.tmp2[ which(outcome == "TRUE")])
plot(mean(var2))

var2 <- readAll(var2)
# readAll() is called to read in all data so that it will be saved in the .rda file Otherwise "raster" will use shortcuts to save memory and results lost when seesion is ended.

save(var2, file=paste(stackpath,"europeStack.rda", sep="/"))
