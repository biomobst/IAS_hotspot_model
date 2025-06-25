# Documentation of hotspot model trial for NIS in 2025

## Summary 

This folder documents the re-run of the original IAS hotspot model with new data for 80 non-indigenous species detected by the genetic monitoring programs ARMS MBON (Pagnier et al 2025) and the Swedish National port monitoring (Sundberg et al 2024). 


The purpose of the models is to match currently known distribution range of marine alien species with its potential areas of suitable habitat in European coastal waters. To this end we analysed 80 species which have been detected outside their native distribution range. From the input data we were able to produce individual models for 52 species, which you can find in this repo. 

## Download data from GBIF and pseudoabsences
Species data are downloaded from the input file “data.table.2025.csv” with the script ”gbif.download.2025.r”. Taxon names and Taxon keys are used to download occurrences (occ_download). A link to the data for download is sent to your email. <br />
<br />
Data filtration takes place in the scriptet ”gbif.download.check.2025.r” , first by excluding erroneous coordinates with the function clean_coordinates() and thereafter by inspecting observations field ”basisOfRecord” and excluding the following classes: "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "MATERIAL_CITATION"  , "MATERIAL_SAMPLE" , "LIVING_SPECIMEN"). The following  classes were included: "HUMAN_OBSERVATION", "OCCURRENCE", "OBSERVATION". <br />
<br />
For each species the cleaned data are saved as.csv file with the following columns: 
"gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "depth", "depthAccuracy","eventDate" <br />


For each species a csv file is generated with the following format 
write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F)

In addition to the files, the script creates plots of finds on a world map, partly for each species individually and partly a plot with all species in different colors to give an image of where in the world finds of invasive species have been made.

The last step is to create pseudoabsences for presence-absence modeling. A plot of all downloaded observations (testplot.cleanput.jpg) shows that observations are unevenly distributed globally and that one can assume that there is a greater chance that someone found a species in an area where other species are reported. 

From the collected downloaded, cleaned, data (first download) 10,000 observations are sampled which are used as pseudoabsences under the assumption that the distribution of all findings in "cleanput" can be considered a reasonable model for "sampling effort".
Of the 10,000 findings, many will, in later stages, turn out to be at points where environmental data is missing, but in total about 1,500 pseudoabsences were useful, which fairly balanced the positive findings.

In the tables "species stats" and "species stats no chlora" there is species-specific information on the number of finds in GBIF that meet the criteria above, the number of finds with unique coordinates and the number of these for which environmental information is available. Data on chlorophyll is more limited and missing in smaller bodies of water, which is why a model based on data with this variable is based on fewer observations and fewer pseudoabsences. 

Information about which files are to be used as input in models for different species is gathered in the file "data.table.csv" which has the following format:

species	present.data	absence.data	pseudoabsence.data
Elodea nuttallii	Elodea nuttallii.csv	Elodea nuttallii.csv	pseudoabsences.csv
Nymphoides peltata	Nymphoides peltata.csv	Nymphoides peltata.csv	pseudoabsences.csv
Pistia stratiotes	Pistia stratiotes.csv	Pistia stratiotes.csv	pseudoabsences.csv

The script that loads data into the model can use different files for positive and negative findings is a inherited function from the previous project. In this project, positive and negative findings are read from the same file and the same pseudoabsences are used for all species, but it is possible, for example, to enter negative findings manually or have different "background models" for different groups of organisms by choosing different pseudoabsence files.

The selection of points in the pseudoabsence model has not been filtered to remove duplicates with the same coordinates before sampling. This means that areas with many duplicates may be overrepresented. Filtering of duplicates was not done until the last run, otherwise it would of course have been done here as well.

The "model" used to handle differences in sampling effort is one of the things that can be further developed in the next project.

Download and formatting of environmental data layers

We downloaded and formatted the following environmental data layers from https://neo.gsfc.nasa.gov/ and generated consistent continental data layers covering land and sea at 5 arcmin resolution for the following variables:

•	Sea surface temperature (SST)
o	Annual mean SST in C
o	Annual max SST in C
o	Annual min SST in C
o	Annual amplitude SST in C
•	Chlorofyll (Chla) 
o	mean Chla of most productive month in mg/m3
o	min Chla of most productive month in mg/m3
o	max Chla of most productive month in mg/m3
o	amplitude Chla of most productive month in mg/m3
•	Sea surface salinity (SSC)
o	Annual mean SSS in PSU

SST- och Chla layers are processed in the script ”Merge NASA GEOTIFF.R”, but also check which file was used to read the salinity values. We used the source https://neo.gsfc.nasa.gov/archive/geotiff/ together with the script ”nasa file download.r” according to the formats:
https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/AQUA_MODIS.YYYYMMDD – YYYYMMDD.L3m.MC.CHL.chlor_a.9km.nc
and 
https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/AQUA_MODIS. YYYYMMDD _ YYYY.MM.DD.L3b.MO.SST.nc


Model test and projection

Modellling was performed using the script invasive.data.masterscript.with CV.2022.r  with help functions library script SEanalytics.functionsNEW.r

The script includes the following steps
•	Filepaths, filenames, an suffixes
•	Read present, absent and pseudoabsence points and extract environmental data
•	Preparing iterations 
•	Estimate weight of predictors with MCMC algortihm
•	Make plots of variables’ weight and correlation with species observations
•	Train random forest models
•	Extract C50 rules 
•	Calculate and plot ROC curves from cross validation 
•	Spatial prediction of species presence and plotting individual maps 
•	Plot map stacks with average probability
•	Plot maps that combine cumulative average probability for all species with traffic layers

Filepaths, filenames, an suffixes
For the most important folders, three variants are created with different "suffixes" for results with and without the chlora variable.

Read present, absent and pseudoabsence points and extract environmental data
The filtered data from GBIF is loaded. Duplicates are then filtered out if they have the same coordinates and occurrence status. Next, environmental data is extracted from the raster stack that has the "correct" suffix. Observations without complete environmental data are filtered out. In the loop, a table is also created with statistics for each species on the number of positive and negative findings, including pseudoabsences. Species with fewer than 5 unique, complete, findings are filtered out.

Preparing iterations 
The analysis is performed as 5 replicates of a 5x5 cross-validation. For each replicate data are permuted. Then, positive and negative findings are sampled, individually, to belong to one of five possible (other values at the CV level are possible) sets, which corresponds to the iteration in the cross-validation when this finding is to constitute the test set.

Since many points are close to each other and would lead to overestimation of predictive power if one allowed nearby points to be included in both the test and training sets, the coordinates are rounded as follows:
  lonmin <-   10*floor(my.data$Lon/10)
  lonmax <-   10*ceiling(my.data$Lon/10)
  latmin <-   10*floor(2*my.data$Lat/10)/2
  latmax <-   10*ceiling(2*my.data$Lat/10)/2
and each observation is given an observation area whose name is given by the above. When you then sample data for training and test sets, observation areas are sampled, not individual observations. For some species, all findings ended up in a few observation areas.

In order to carry out a 5x cross-validation, there must be findings in at least 5 different areas, otherwise the algorithm crashes. Since a meaningful estimate of predictive power, and also confidence in the maps, is doubtful if the positive findings come from fewer than five areas, such species are filtered out.

The size of the rasters and the very principle of grouping findings that may be considered too close to each other for meaningful ROC analysis can be debated.

Information about the iterations is saved in a .rda file and is then used partly in the variable selection algorithm and partly when the Random Forest model is trained for cross-validation.

Estimate weight of predictors with MCMC algorithm
The weight of the variables is estimated with MCMC feature selection according to the same principle as in earlier model trials by Bergkvist et al (2020). Since the number of variables is small, however, we do not use the possibility of making a Random Forest model with only significant variables, but MCMC is only used to get a measure of the usefulness of the variables. The MCMC algorithm is run 25 times with different parts of the dataset to get a feel for how sensitive the weight of the variables (RI index) is to the selection of data. However, as the method is implemented, nearby points may be included in the same dataset and the weight of the variables may be affected by overtraining. This should be handled in the next project, for example by filtering out nearby finds.

The MCMC algorithm is time-consuming and the 25 iterations are carried out on different nodes with the program package parallel{}.

Make plots of variables’ weight and correlation with species observations 
In part, a plot is made per species that shows the weight of the variables (RI index). A couple of random variables have been added to the model to prevent it from crashing. In the plot, a horizontal line has been inserted corresponding to the RI index that the model considers to be "better than chance". In addition, plots are made for each variable where x is the measured value of the predictor variable divided into 10 equally sized "bins" and the axis is the proportion of positive findings given that x lies in this bin.

Train random forest models
A model is trained for each sub-dataset according to the file created under “prepare iterations”. For these models, only the "probability that the observation has the status present" is saved for each observation. In addition, a model is trained with the entire input data that is saved and used for spatial prediction.

In this implementation, Random Forests are run sequentially. The time gain of running in parallel on different nodes has not been that great, as the number of variables is small, and it has not been worth the time to modify the script for parallel training of models. However, this could be done easily in the future. 


Extract C50 rules 
Test of a simple algorithm to create simple "rules" that describe the relationships between predictor variables and occurrence status. The output of this is saved as a text file for each species. The algorithm can be seen as an attempt to explain what happens inside the "decision trees" even if it is a different decision tree algorithm. The function has been tested without any optimization whatsoever, but it might be interesting to compare with the plots made above and as a "demo" of a track to follow up in future projects.

Calculate and plot ROC curves from cross validation 
ROC kurvor beräknas baserat från Random Forest förningarna ovan. Kurvorna baseras på den sannolikhet som predikerats för varje observation när den utgjorde testdata, och alltså med en modell som inte bygger på punkter i denna rektangel. Eftersom Random Forest algoritmen körts i 5 repetitioner av en 5x korsvalidering får 5 ROC kurvor med olika AUC. Alla fem kurvorna visas som en linje i plotten och medelkurvan som en grön polygon. Hur stor skillnad det är mellan AUC värdena från olika repetitioner ger en indikation på hur känslig modellen är för data.

Spatial prediction of species presence and plotting individual maps (Sweden and global)
ROC curves are calculated based on the Random Forest simulations above. The curves are based on the probability predicted for each observation when it constituted the test data, and thus with a model not based on points in this rectangle. Since the Random Forest algorithm was run in 5 repetitions of a 5x cross-validation, 5 ROC curves with different AUCs are obtained. All five curves are shown as a line in the plot and the average curve as a green polygon. How much difference there is between the AUC values from different repetitions gives an indication of how sensitive the model is to the data.

The Random Forest model trained with the entire dataset is used for prediction with a raster stack corresponding to the data used to train the model. The results are first saved as .rda files and, in a later step, as GeoTIFF. There is room here to make the script more uniform and skip the first step.

Plot map stacks with average probability
For each experiment, a raster stack is made with predicted probability from all species for which the model succeeded in creating such a map. An average value is calculated with the mean() function. There are also maps with logarithmic probability and weighted logarithmic probability. This is a legacy from the previous project by Bergkvist et al (2020). 

Plot maps that combine cumulative average probability for all species with traffic layers
Raster layers with traffic data from 2016 that were used in the previous project (Bergkvist et al 2020) are read in and then the raster map corresponding to average probability is cropped and resampled to the same extent and resolution as the traffic map. The area being analyzed is given by the extent of the map used by Bergkvist et al (2020). 

A composite plot is made with basically the same code as in the previous project. The plot has four panels with


•	Average probability for all species given type of model (with or without chlorophyll and without chlorophyll but only points where chlorophyll data is available)
•	Traffic data from 2016. The difference to the previous analysis is that the ceiling for traffic intensity was set to 1000 instead of 10000 (values above the ceiling are set to the ceiling). The reason for truncating high values is that otherwise you don't see the traffic lanes in the Gulf of Bothnia, which are much less intensive than the major waterways in the southern Baltic Sea.
•	A map where the traffic data is superimposed on the model results with a different color.
•	A map where traffic data and model results have been added. During the addition, a scale parameter is used which is set arbitrarily so that the traffic lane will appear.


 
### References

Daraghmeh N, Exter K, Pagnier J, et al (2024) A long-term ecological research data set from the marine genetic monitoring programme ARMS-MBON 2018-2020. Molecular Ecology Resources. doi: https://doi.org/10.1111/1755-0998.14073

Pagnier J, Daraghmeh N, Obst M (2024) Using the long-term genetic monitoring network ARMS-MBON to detect marine non-indigenous species along the European coasts. Biological Invasions 27, 77 (2025). https://doi.org/10.1007/s10530-024-03503-2

Sundberg P, Obst M, Panova M (2024) DNA-baserad övervakning av arter i akvatisk miljö - verifiering och tillämpning. Naturvårdsverket rapport. pp 44. Rapport 7157. ISBN 978-91-620-7157-8
