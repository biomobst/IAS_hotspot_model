# open files downloaded from https://neo.gsfc.nasa.gov/archive/geotiff/
require(raster)
require(rgdal)
library(sf)
require(fBasics)

path <- "~/path/to/project/Data2022/NASA"
# set one path per class of downloaded data.
rasterpaths <- paste(path, c("CHLORA tiff/", "SST tiff/"), sep="/")

processpath <-paste(path, "processlayers" ,sep="/")
layersclass <- c("CHLORA", "SST")
#Chech that files are there
for(mypath in rasterpaths){
lista.ras<- Sys.glob(paste(mypath,"*.TIFF",sep=""))
print( lista.ras)
}

#get data for same month all years
months <- c("01","02","03","04","05", "06","07","08","09","10","11","12")
my.maps <- list()
#for(mypath in rasterpaths){
for(i in 1:2){
  mypath <- rasterpaths[i]
  my.maps[[layersclass[i]]] <- list()
  
  for(month in months){
    lista.ras<- Sys.glob(paste(mypath,"*",month,".TIFF",sep=""))
    print( lista.ras)
     mystack <- stack(lista.ras)

     tf <- writeRaster(mystack, paste(processpath,"/",layersclass[i],
                                "_",month,".TIFF", sep=""),format="GTiff",overwrite=TRUE)
  }
}
gc()
for(i in 1:2){
  mypath <- rasterpaths[i]
  for(month in months){
    gc()
    filename <- paste(processpath,"/",layersclass[i],"_",month,".tif", sep="")
    print(filename)
    filename2 <- paste(processpath,"/",layersclass[i],"_",month,".mod.tif", sep="")
        mystack <- stack(filename)
        # define that extreme values are NAs
    mystack[mystack>249]<- NA
    writeRaster(mystack, filename2, format="GTiff",overwrite=TRUE)
    print(filename2)

  }
}

# Options to plot maps

#for(names(mystack))
#map <- raster(lista.ras[[2]])
#par(mfrow =c(1,1))
#plot(map)

#n <- length(names(mystack))
# names(my.maps)
#par(mar = c(0,0,0,5))
#par(mfrow =c(4,3))
#for (i in 1:12){
#  plot(my.maps[[2]][[9]])
#}
#par(mfrow =c(4,3))
#for (i in 1:12){
#  plot(my.maps[[1]][[i]])
#}
#par(mfrow =c(1,1))
#plot(my.maps[[2]][[7]])


#par(mfrow =c(2,2))
#plot(my.maps[[1]][[7]], xlim = c(-30,70), ylim = c(20,70))
#title(main = names(my.maps[1]))
#plot(my.maps[[1]][[7]], xlim = c(-150,-50), ylim = c(20,70))
#title(main = names(my.maps[1]))


#plot(my.maps[[2]][[7]], xlim = c(-30,70), ylim = c(20,70))
#title(main = names(my.maps[2]))
#plot(my.maps[[2]][[7]], xlim = c(-150,-50), ylim = c(20,70))
#title(main = names(my.maps[2]))


#par(mfrow =c(2,1))
#plot(my.maps[[1]][[7]], xlim = c(-10,0), ylim = c(36,44))
#title(main = names(my.maps[1]))
#plot(my.maps[[2]][[7]], xlim = c(-10,0), ylim = c(36,44))
#title(main = names(my.maps[2]))

################################################
# Prepare average of same month from all years. Disregard NAs
my.meansmaps <- list()
my.meansmaps[["CHLORA"]] <- list()
my.meansmaps[["CHLORA"]][["max"]] 
for(i in 1:2){
  filenames <- paste(processpath,"/",layersclass[i],"*",".mod.tif", sep="")
  
   lista.ras<- Sys.glob(filenames)
   mystack <-stack(lista.ras)
  
  map <- mean(mystack, na.rm = T)
  filename<- paste(processpath,"/",layersclass[i],"_",".MEAN.tif", sep="")
  writeRaster(map, filename, format="GTiff",overwrite=TRUE)
  
  map <- max(mystack, na.rm = T)
  filename<- paste(processpath,"/",layersclass[i],"_",".MAX.tif", sep="")
  writeRaster(map, filename, format="GTiff",overwrite=TRUE)
  
  map2 <- min(mystack, na.rm = T)
  filename<- paste(processpath,"/",layersclass[i],"_",".MIN.tif", sep="")
  writeRaster(map2, filename, format="GTiff",overwrite=TRUE)
  
  map3 <- map-map2
  filename<- paste(processpath,"/",layersclass[i],"_",".AMP.tif", sep="")
  writeRaster(map3, filename, format="GTiff",overwrite=TRUE)

}
  
chlorastack <-stack(my.maps[[1]])
chloramean<- mean(chlorastack, na.rm = T)
plot(chloramean)

tempstack <-stack(my.maps[[2]])
mystack[mystack>249]<- NA

tempmean<- mean(tempstack, na.rm = T)
plot(tempmean)
tempmax<- max(tempstack, na.rm = T)
tempmin<- min(tempstack, na.rm = T)
tempamp <- tempmax - tempmin

par(mfrow =c(2,2))
par(mar = c(0,0,4,4))
colorsBlue <- seqPalette(n= 20, name = c( "Blues") ) 
#heatColors <- rampPalette(n = 10, name = c("blue2red")) 
heatColors <- rev(seqPalette(20, name = c("YlGnBu")))
colorsGreen <- seqPalette(n=20, name = c("Greens"))
#ramp

for (i in 1:2){
  filename<- paste(processpath,"/",layersclass[i],"_",".MEAN.tif", sep="")
  plot(raster(filename), ylim = c(-80, 80), col = list(colorsGreen, heatColors)[[i]])
  title(main = layersclass[i])
  
  filename<- paste(processpath,"/",layersclass[i],"_",".MIN.tif", sep="")
  plot(raster(filename), ylim = c(-80, 80),col = list(colorsGreen, heatColors)[[i]]) 
  title(main = "Min")
  
  filename<- paste(processpath,"/",layersclass[i],"_",".MAX.tif", sep="")
  plot(raster(filename), ylim = c(-80, 80), col = list(colorsGreen, heatColors)[[i]] )
  title(main = "Max")
       
  filename<- paste(processpath,"/",layersclass[i],"_",".AMP.tif", sep="")
  plot(raster(filename), ylim = c(-80, 80), col = colorsBlue)
  title(main = "Amp")
  
}

par(mfrow = c(1,1))
filename<- paste(processpath,"/",layersclass[1],"_",".MAX.tif", sep="")
plot(raster(filename), ylim = c(-80, 80), col = list(colorsGreen, heatColors)[[1]] )
title(main = paste( layersclass[1],"Max"))

par(mfrow = c(1,1))
filename<- paste(processpath,"/",layersclass[1],"_",".MAX.tif", sep="")
plot(raster(filename), ylim = c(30, 50),xlim =c(-10, 60) ,col = list(colorsGreen, heatColors)[[1]] )
title(main = paste( layersclass[1],"Max"))

par(mfrow = c(1,1))
filename<- paste(processpath,"/",layersclass[2],"_",".MAX.tif", sep="")
plot(raster(filename), ylim = c(30, 50),xlim =c(-10, 60) ,col = list(colorsGreen, heatColors)[[2]] )
title(main = paste( layersclass[2],"Max"))

par(mfrow =c(4,3))
for (i in 1:12){
  plot(chlorastack[[i]])
}


########################################################
#Add a layer with salinity (SSS)
##############################################################
i <- 2
# select one of the NASA datalayers as mask
filename<- paste(processpath,"/",layersclass[i],"_",".AMP.tif", sep="")
mask <- raster(filename)
extent(mask)
plot(mask)

# Load the file with rasterstack from last project and select salinity layer
load("~/path/to/project/Data2022/Rasterstacks/globalStack.rda")
names(var)

#[80] "Present.Surface.Salinity.Max"                         
#[81] "Present.Surface.Salinity.Mean"                        
#[82] "Present.Surface.Salinity.Min"                         
#[83] "Present.Surface.Salinity.Range" 

mean.salinity <- var[[81]] 
plot(mean.salinity)
#Define that Freshwater salinity: < 0,05 %
extent(mean.salinity)
dim(mask)
dim(mean.salinity)
mean.salinity <- resample(mean.salinity, mask, method="bilinear")
dim(mean.salinity)
mean.salinity[is.na(mean.salinity)]<- 0.5 # Set a defalult salinity where there is no data. (also on land)
plot(mean.salinity)
mean.salinity[is.na(mask)]<- NA # Define that there is no salinity where there is no SST data. Assumed to be land.
plot(mean.salinity)
tf <- writeRaster(mean.salinity, paste(processpath,"/mean.salinity.TIFF", sep=""),format="GTiff",overwrite=TRUE)

# Prepare a new rasterstack that also includes salinity
layers <- c(paste(processpath,"/",layersclass[2],"_",".AMP.tif", sep=""),
paste(processpath,"/",layersclass[2],"_",".MIN.tif", sep=""),
paste(processpath,"/",layersclass[2],"_",".MAX.tif", sep=""),
paste(processpath,"/",layersclass[2],"_",".MEAN.tif", sep=""),
paste(processpath,"/",layersclass[1],"_",".MAX.tif", sep=""),
paste(processpath,"/mean.salinity.TIFF", sep =""))

rasterstack.global.2022 <- stack(layers)
names(rasterstack.global.2022)[6]<- "SALINITY"
tf <- writeRaster(rasterstack.global.2022, paste(processpath,"/rasterstack.global.2022.TIFF", sep=""),format="GTiff",overwrite=TRUE)

e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax 
rasterstack.Europe.2022 <- crop(rasterstack.global.2022, e)	
plot(rasterstack.Europe.2022)
names(rasterstack.Europe.2022)#[6]<- "SALINITY"

tf <- writeRaster(rasterstack.Europe.2022, paste(processpath,"/rasterstack.Europe.2022.TIFF", sep=""),format="GTiff",overwrite=TRUE)

# store the layernames as a separate file, since the information is lost when saving rasterstack
layernames <- names(rasterstack.Europe.2022)
save(layernames, file = paste(processpath,"/layernames.rda" ,sep =""))


#############################################################################
### make stack without chlora but only sites where chlora is available

tf <- stack(paste(processpath,"/rasterstack.Europe.2022.TIFF", sep=""))

load(file = paste(processpath,"/layernames.rda" ,sep =""))
names(tf) <- layernames

plot(tf[[5]])
tf[[6]][is.na(tf[[5]])]<- NA
tf[[1]][is.na(tf[[5]])]<- NA
tf[[2]][is.na(tf[[5]])]<- NA
tf[[3]][is.na(tf[[5]])]<- NA
tf[[4]][is.na(tf[[5]])]<- NA

plot(tf[[6]])
tf <- tf[[-5]]
layernames <- names(tf)
writeRaster(tf, paste(processpath,"/rasterstack.Europe.2022.no.chlora.delim.TIFF", sep=""),format="GTiff",overwrite=TRUE)
save(layernames, file = paste(processpath,"/layernames.no.chlora.delim.rda" ,sep =""))


#### same as above at global scale
tf <- stack(paste(processpath,"/rasterstack.global.2022.TIFF", sep=""))

load(file = paste(processpath,"/layernames.rda" ,sep =""))
names(tf) <- layernames

plot(tf[[5]])
tf[[6]][is.na(tf[[5]])]<- NA
tf[[1]][is.na(tf[[5]])]<- NA
tf[[2]][is.na(tf[[5]])]<- NA
tf[[3]][is.na(tf[[5]])]<- NA
tf[[4]][is.na(tf[[5]])]<- NA

plot(tf[[6]])
tf <- tf[[-5]]
layernames <- names(tf)
writeRaster(tf, paste(processpath,"/rasterstack.global.2022.no.chlora.delim.TIFF", sep=""),format="GTiff",overwrite=TRUE)
save(layernames, file = paste(processpath,"/layernames.no.chlora.delim.rda" ,sep =""))



#################################################################################
###################################################################################
#make stack without chlora
tf1 <- stack(paste(processpath,"/rasterstack.global.2022.TIFF", sep=""))
tf2 <- stack(paste(processpath,"/rasterstack.Europe.2022.TIFF", sep=""))
load(paste(processpath,"/layernames.rda" ,sep =""))
layernames
names(tf1) <- layernames
names(tf2) <- layernames

tf1 <- tf1[[ c("SST_.AMP","SST_.MIN","SST_.MAX","SST_.MEAN", "SALINITY" )]]
tf2 <- tf2[[ c("SST_.AMP","SST_.MIN","SST_.MAX","SST_.MEAN", "SALINITY" )]]

writeRaster(tf1, paste(processpath,"/rasterstack.global.2022.no.chlora.TIFF", sep=""),format="GTiff",overwrite=TRUE)
writeRaster(tf2, paste(processpath,"/rasterstack.Europe.2022.no.chlora.TIFF", sep=""),format="GTiff",overwrite=TRUE)
layernames <- names(tf1)
save(layernames, file = paste(processpath,"/layernames.no.chlora.rda" ,sep =""))



########