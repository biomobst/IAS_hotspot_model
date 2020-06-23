This folder contains all individual SDM results (species level) in raster format (GeoTIFF).

Habitat suitability values are expressed with two types scales (see below for explanation). 

All maps are provided as European projections.

## Scaling ##

linear.prob. = The probability of class “present”. Obtained with function predict(rasterstack, RandomForestmodel, type="prob")
 
log.prob. = Log10(linear.prob + 0.001) # to avoid division by zero

Range linear.prob. [0 -> 1]
Range log.prob. [-3 -> 0]
