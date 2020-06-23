This folder contains all summarised SDM results in raster format (GeoTIFF). The maps show suitable habitat for the below ecological groupings of species and with 4 different levels of weighting implemented (based on individual risk assessments). 

Habitat suitability values are expressed with two types of logarithmic scales (see below for explanation). 

All maps are provided as European projections.

## Species groupings

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


## Weighting  and cutoffs

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

## Scaling

For each species a matrix is calculated as 3+(log.prob).  This takes value 3 when probability is 1 (thus log10(prob = 1) and 0 when probability is 10exp-3  (or lower).

The summary matrix is the sum of the matrix  for all species in the group given by the filename.

If filename contains “raw” the sum is unweighted.

raw.log = Sum(  for species in group.of.species) { 3+log.prob[species]  }  )
if filename contains “weighted” the summary matrix is weighted by “samlat riskutfall”

weighted.log =  Sum(  for species in group.of.species) { (3+log.prob[species]   )* Samlat riskutfall[species]  }  )
