library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
#library(rgdal)

library(raster)
library(ncdf4)
library(RNetCDF)

projectpath <- "~/Dokument/Projekt/HAV2025"
dir <- "~/Dokument/Projekt/HAV2025/data/Biooracle.download/"
dataset_scenarios <- c( "baseline" ,"ssp119" ,  "ssp126"  , "ssp245"   ,"ssp370" ,  "ssp460"  , "ssp585"  )

#for(sel.sen in c(2:7)){
  sel.sen <- 1
  scenario <- dataset_scenarios[[sel.sen]]
  #  outdir <- paste(dir,"/datalayer.tiff/",scenario,"/",sep="")
  
stackpath <- paste(dir,"rasterstacks/",scenario,"/",sep="")
rasterpath <- paste(dir,"datalayer.nc/",scenario,sep="")

#paste(path, c("data/Biooracle_datalager/datalager.tiff"), sep="/")

lista.ras<- Sys.glob(paste(rasterpath,"/*.nc",sep=""))
layers <- rast(lista.ras)  

# Use the first layer from your NetCDF stack as the template:
#layers()
#mystacknc <- stack(lista.ras[1])
#mystacknc <- ncdf4::nc_open(lista.ras[1])

template <- layers[[1]]

# Load and reproject the land polygons to the template CRS:
land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, crs(template))

# Convert to a SpatVector:
land_vect <- vect(land)
plot(land_vect)

plot(template)
# Rasterize the land polygons onto the template.
# Use background = NA so that cells not covered by a polygon remain NA.
land_mask <- rasterize(land_vect, template, field = 1, background = NA)

# Now force any cell that got a value (i.e. where land is present) to 1:
land_mask[!is.na(land_mask)] <- 1

# Check the unique values:
print(unique(values(land_mask)))
# Expect something like: [1]  1  NA

masked_list <- lapply(1:nlyr(layers), function(i) {
  this_layer <- layers[[i]]
  mask(this_layer, land_mask, maskvalue = 1)
})
masked_layers <- rast(masked_list)

filled_list <- lapply(1:nlyr(masked_layers), function(i) {
  this_layer <- masked_layers[[i]]
  focal(this_layer, w = 3, fun = function(x) {
    if (is.na(x[5])) {
      mean(x, na.rm = TRUE)
    } else {
      x[5]
    }
  })
})
filled_layers <- rast(filled_list)

# PLOTS TO CHECK RESULTS

plot(land_mask, main = "Binary Land Mask (1 = Land, NA = Water)")
# Overlay the original land polygons for reference
plot(land_vect, add = TRUE, border = "red")

plot(masked_layers[[1]], main = "Masked Layer (Land = NA)")
# Overlay the land polygons
plot(land_vect, add = TRUE, border = "blue")

plot(filled_layers[[1]], main = "Filled Layer (Coastal Interpolated)")
# Overlay the land polygons to see the boundary
plot(land_vect, add = TRUE, border = "green")

par(mfrow = c(1, 3))
plot(layers[[1]], main = "Original Layer")
plot(masked_layers[[1]], main = "Masked Layer")
plot(filled_layers[[1]], main = "Filled Layer")
par(mfrow = c(1, 1))


summary(values(template))
summary(values(filled_layers[[1]]))
writeRaster(filled_layers, paste(stackpath ,"filled_layers_new.tif",sep=""), overwrite = TRUE, filetype = "GTiff")
saveRDS(filled_layers, paste(stackpath ,"filled_layers.rds",sep=""))

########## load and create stacks for use
#filled_layers <- readRDS(paste(stackpath ,"filled_layers_new.tif",sep=""))
mystack <- stack(paste(stackpath ,"filled_layers_new.tif",sep=""))
plot(mystack)
e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax 
rasterstack.filled.layers.Europe.2025 <- crop(mystack, e)
plot(rasterstack.filled.layers.Europe.2025)

filename <- paste(stackpath,"/","Biooracle.filled.layers.global2025.tif", sep="")
writeRaster(mystack, filename, format="GTiff",overwrite=TRUE)

filename2 <- paste(stackpath,"/","Biooracle.filled.layers.Europe2025.tif", sep="")
writeRaster(rasterstack.filled.layers.Europe.2025, filename2, format="GTiff",overwrite=TRUE)


layernames <- names(mystack)
save(layernames, file = paste(stackpath,"/layernames.filled.rda" ,sep =""))
# get full layernames
load(paste(stackpath,"/layernames.rda" ,sep =""))

names(mystack) <- layernames