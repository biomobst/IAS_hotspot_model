require(raster)
#require(rgdal)
library(sf)
require(fBasics)
        
projectpath <- "~/Dokument/Projekt/HAV2025"
dir <- "/home/gunnarandersson/Dokument/Projekt/HAV2025/data/Biooracle.download/"
dataset.scenarios <- c( "baseline" ,"ssp119" ,  "ssp126"  , "ssp245"   ,"ssp370" ,  "ssp460"  , "ssp585"  )

for(sel.sen in c(1:7)){
  #sel.sen <- 1
  scenario <- dataset_scenarios[[sel.sen]]
#  outdir <- paste(dir,"/datalayer.tiff/",scenario,"/",sep="")
  
rasterpath <- paste(dir,"/datalayer.tiff/",scenario,"/",sep="")
stackpath <- paste(dir,"/rasterstacks/",scenario,"/",sep="")
if (!dir.exists(stackpath)) dir.create(stackpath) 

lista.ras<- Sys.glob(paste(rasterpath,"/*",".tif",sep=""))
mystack <- stack(lista.ras)
for(i in lista.ras){
  print(i)
  #plot(raster(i))
  gc()
}
filename <- paste(stackpath,"/","Biooracle.global2025.tif", sep="")
writeRaster(mystack, filename, format="GTiff",overwrite=TRUE)

e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax 
rasterstack.Europe.2025 <- crop(mystack, e)
plot(rasterstack.Europe.2025)

filename2 <- paste(stackpath,"/","Biooracle.Europe2025.tif", sep="")
writeRaster(rasterstack.Europe.2025, filename2, format="GTiff",overwrite=TRUE)


layernames <- names(mystack)
save(layernames, file = paste(stackpath,"/layernames.rda" ,sep =""))

}