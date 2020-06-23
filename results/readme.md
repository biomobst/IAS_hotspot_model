# DOCUMENTATION, version 2020-05-24
 
## TargetSpeciesList_20200107

This table contains the information on the species included in the study at this stage as well as associated ecological information and risk groupings, etc.

## FOLDER: Maps SPECIES 

This folder contains all individual SDM results (species level) as PNG files. The maps showi suitable habitat for each species as European projection (name.Eur.png) and Swedish projection (name.Swe.png).

The scale is linear from white to dark brown, with 0 indicating lowest probability to find suitable habitat, and 1 indicating highest  probability to find suitable habitat. 

Training points (based on known observations) are indicated with red circle in the maps.

## FOLDER: Raster files SPECIES 

This folder contains all individual SDM results (species level) in raster format (GeoTIFF).

Habitat suitability values are expressed with two types scales (see below for explanation). 

All maps are provided as European projections.

### Scaling

linear.prob. = The probability of class “present”. Obtained with function predict(rasterstack, RandomForestmodel, type="prob")
 
log.prob. = Log10(linear.prob + 0.001) # to avoid division by zero

Range linear.prob. [0 -> 1]
Range log.prob. [-3 -> 0]

## FOLDER: Raster files GROUPS 

This folder contains all summarised SDM results in raster format (GeoTIFF). The maps show suitable habitat for the below ecological groupings of species and with 4 different levels of weighting implemented (based on individual risk assessments). 

Habitat suitability values are expressed with two types of logarithmic scales (see below for explanation). 

All maps are provided as European projections.

### Species groupings

- all species (all species in the masterfile included)
- Marine (only marine species included, i.e. we omitted freshwater species that usually migrate through rivers and lakes across the European continent. Often with origin in the Ponto-Caspain region)
- Freshwater (only freshwater tolerant species usually migrating through rivers and lakes across the European continent. Often with origin in the Ponto-Caspain region)
- Plankton (only species that are as adults in zooplankton and phytoplankton)
- Benthos (only species that are as adults in zoobenthos and phytobenthos)
- Zooplankton and zoobenthos (only species that are as adults in zooplankton and zoobenthos)
- Phytoplankton and phytobenthos (species that are as adults phytoplankton and phytobenthos)
- Phytobenthos (only phytobenthos species)
- Phytoplankton (only phytoplankton species)
- Zooplankton (only zooplankton species)
- Zoobenthos (only zoobenthos species)

### Weighting  and cutoffs

We calculated weighted and unweighted sums when when summarising the individual projections into group projections. The purpose is to give species different impacts on the grouped maps based on HaV's risk assessment.   

We also applied different cut-offs and calculated weighted and unweighted sums of species where predicted probability of presence exceeded a threshold of 0.1, 0.33 and 0.5.

“sum species” in filename.

Cutoff given by filename

“raw” in filename = unweighted. 
“weighted” in filename = risk-weighted by number “samlat riskutfall” in the Overview file.

Weighted sum=:  Sum(  for species in group.of.species) { If(linear.prob[species] >= cutoff){1}else{0} ) * Samlat riskutfall[species] ) }  )

raw sum=:  Sum(  for species in group.of.species) {  If(linear.prob[species] >= cutoff){1}else{0} )   }  )

group of species defined by environment or salinity group

EXAMPLE: sum species weighted 0.5 marine environment:
Sum of species in the group "marine environment". Raster values only included if higher than 0.5 and then multiplied by risk factor for that species.

Range "sum species

### Scaling

For each species a matrix is calculated as 3+(log.prob).  This takes value 3 when probability is 1 (thus log10(prob = 1) and 0 when probability is 10exp-3  (or lower).

The summary matrix is the sum of the matrix  for all species in the group given by the filename.

If filename contains “raw” the sum is unweighted.

raw.log = Sum(  for species in group.of.species) { 3+log.prob[species]  }  )
if filename contains “weighted” the summary matrix is weighted by “samlat riskutfall”

weighted.log =  Sum(  for species in group.of.species) { (3+log.prob[species]   )* Samlat riskutfall[species]  }  )

## FOLDER: Training points 

This folder contains the Pseudo-absence points used in the model as well as the species presence and absence data as csv files.

## FOLDER: Confusion matrices 

This folder contains the confusion matrices or out-of-bag cases that were used as estimates of predictive performance on new cases.

## ALGORITHM 

For each species:

Read file(s) with observed presence and absence points

Read file with pseudo absence points

Combine tables: make “presence”/”absence” a factorial variable

Obtain environmental data for each point

Select training set

Calculate class weights for classes present and absent

  w.pres = length(observations(which(class == "present")))/length(observations)

  w.abs = 1

train random forest model with training data and given classweights

Calculate OUB error.

save model + performance parameters

As the performance of the models were calculated on out-of-bag cases (OUB error) upon training, we extracted this from the stored models and printed in attached file.

