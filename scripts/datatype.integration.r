source("SEanalytics.functions.r")

require(raster)
require(rgdal)
require(zip)

#######

#rastermappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters"
rastermappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters/predict rasters apr 2020"
trafficpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/traffic layers/2016_All_shiptypes_AIS_Shipping_Density"
#Plotpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/plots/combiplots" 
Plotpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/plots/plots apr 2020"
groups <- c("*sum species*","*weighted.log*","*raw.log*")

#group <- "*raw.log.map.all species*"
#group <- "*mean of log10 prob all.tif"

#
lista.ras <- c("C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/predict.rasters/log.prob.Anteaeolidiella lurana.tif")
#group <- "*raw.log.map.freshwater.tif"
for(group in groups){
  #brk <- c(seq(0, 1,by=0.1),1.1)
  lista.ras <- Sys.glob(paste(rastermappath,group,sep="/"))
  plotvar <- stack(lista.ras)
  print(summary(plotvar[[1]]))
}

lista.ras2 <- Sys.glob(paste(trafficpath,"*.tif",sep="/"))
shipvar <- stack(lista.ras2)
par(mar=rep(1,4))
plot(plotvar[[1]])
plot(shipvar[[1]])
proj4string(shipvar)
#2016 All shiptypes AIS Shipping Density1

shape3 <- readOGR(dsn=trafficpath,layer="2016 All shiptypes AIS Shipping Density1.tif.vat")
proj4string(shape2)#GRS80 EPSG:4019
shape2 <- spTransform(shape2, CRS("+init=epsg:4326"))
plot(shape2)# just to check that it appears right

new.plotvar <- resample(plotvar[[1]],shipvar[[1]] , method="bilinear")
extent(plotvar[[1]])
extent(shipvar[[1]])
#shipvar2 <- spTransform(shipvar[[1]], CRS("+init=epsg:4326"))
shipvar2 <- projectRaster(shipvar[[1]],crs=CRS("+init=epsg:4326"))
extent(shipvar2)
extent(plotvar[[1]])

new.plotvar <- resample(plotvar[[1]],shipvar2 , method="bilinear")
#colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"#D3D3D3","#D3D3D3","#D3D3D3" )
#colorsorange <- c(seqPalette(n=,12, name = c( "Oranges") ),"#D3D3D3","#D3D3D3","#D3D3D3" )
#colorsdiv <- c(divPalette(n=,12, name = c( "Oranges") ),"#D3D3D3","#D3D3D3","#D3D3D3" )

colorsBlue <- seqPalette(n=,12, name = c( "Blues") )
colorsorange <- seqPalette(n=,12, name = c( "Oranges") ) 
colorsdiv <- divPalette(n=,12, name = c( "BrBG") )

newfile <- paste(Plotpath,"/demoplot.freshwater.png",sep="")
png(newfile,  width = 180, height = 180, units = "mm", res=1200)


plot(new.plotvar, col= colorsBlue)
#arg <- list(at=c(0,0), labels=c("NA","NA"))
plot(shipvar2, add=T, col=colorsorange,alpha=0.5, legend=FALSE)
plot(shape2, add=T)
dev.off()
#################
#plot(shipvar2)
shipvar3 <- shipvar2
shipvar3[is.na(shipvar3)] <- 0

shipvar3[is.na(new.plotvar)] <- NA

shipvar4 <- shipvar3
shipvar4[is.na(shipvar3)] <- 0
shipvar4[shipvar3>10000] <- 10000

combivar <- new.plotvar + shipvar3/5000
combivar[combivar > 3 ]<- 3

newfile <- paste(Plotpath,"/demoplot2.Anteaeolidiella lurana.png",sep="")
png(newfile,  width = 180, height = 180, units = "mm", res=1200)
plot(combivar, col= colorsdiv)
plot(shape2, add=T)
dev.off()

#############################################
#############################################
par(mfrow=c(2,2))
par(mar=c(2,0,2,0))
newfile <- paste(Plotpath,"/demoplot3.freshwater.png",sep="")
png(newfile,  width = 360, height = 360, units = "mm", res=1200)
par(mfrow=c(2,2))
par(mar=c(0,2,6,8))

plot(new.plotvar, col= colorsBlue)
title(main="weighted risk")
plot(shape2, add=T, lwd=0.1, col="lightgrey")

plot(shipvar3, col= colorsorange)
title(main="traffic")
plot(shape2, add=T, lwd=0.1, col="lightgrey")

new.plotvar2 <- new.plotvar
new.plotvar2[new.plotvar2>1.5]<-1.5
plot(new.plotvar2, col= colorsBlue)
#arg <- list(at=c(0,0), labels=c("NA","NA"))
plot(shipvar4, add=T, col=colorsorange,alpha=0.5, legend=FALSE)
title(main="superimposed")
plot(shape2, add=T, lwd=0.1, col="lightgrey")


plot(combivar, col= colorsBlue)
title(main="added by function")
plot(shape2, add=T, lwd=0.1, col="lightgrey")
dev.off()
