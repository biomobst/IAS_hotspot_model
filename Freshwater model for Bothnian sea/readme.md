# Documentation of hotspot model for the Gulf of Bothnia

## Summary 

This folder documents the a customisation of the original IAS hotspot model with the purpose to predict suitable habitat in the northern Baltic Sea for invasive species originating from freshwater systems, i.e. lakes and rivers. The model can be used to estimate the risk of limnic species spreading into or via the Gulf of Bothnia. We analyzed 34 species from different taxonomic groups and were able to produce individual models for 26 species, which you can find in this repo. The maps show that for many species found in northern European freshwater systems there are also suitable habitats in the Gulf of Bothnia. There is thus a risk for freshwater species to spread into the Gulf of Bothnia and other parts of the northern Baltic Sea.  The analysis of the models' predictions in relation to waterways that can spread the species between lakes and the sea and between different parts of the northern Baltic Sea justify monitoring port areas with large estuaries in the Gulf of Bothnia.

## Download data from GBIF and pseudoabsences
Species data are downloaded from the input file “valda.arter.4.csv” with the script ”load.from.gbif.r”. Taxon names and Taxon keys are used to download occurrences (occ_download). A link to the data for download is sent to your email. <br />
<br />
Data filtration takes place in the scriptet ”gbif.download.check.r” , first by excluding erroneous coordinates with the function clean_coordinates() and thereafter by inspecting observations field ”basisOfRecord” and excluding the following classes: "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "MATERIAL_CITATION"  , "MATERIAL_SAMPLE" , "LIVING_SPECIMEN"). The following  classes were included: "HUMAN_OBSERVATION", "OCCURRENCE", "OBSERVATION". <br />
<br />
For each species the cleaned data are saved as.csv file with the following columns: 
"gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "depth", "depthAccuracy","eventDate" <br />
<br />
For each species a csv file is generated with the following format 
write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F) <br />
<br />
In addition to the files, the script creates plots of finds on a world map, partly for each species individually and partly a plot with all species in different colors to give an image of where in the world finds of invasive species have been made. <br />
<br />
The last step is to create pseudoabsences for presence-absence modeling. A plot of all downloaded observations (testplot.cleanput.jpg) shows that observations are unevenly distributed globally and that one can assume that there is a greater chance that someone found a species in an area where other species are reported. <br />
<br />
From the downloaded data files, cleaned data for 10,000 observations are sampled and used as pseudoabsences under the assumption that the distribution of all findings in "cleanput" can be considered a reasonable model for "sampling effort". Of the 10,000 findings, many will, in later stages, turn out to be at points where environmental data is missing, but in total about 1,500 pseudoabsences were useful, which fairly well balanced the positive findings for most species.<br />
<br />
In the tables "species stats" and "species stats no chlora" there is species-specific information on the number of finds in GBIF that meet the criteria above, the number of finds with unique coordinates and the number of these for which environmental information is available. Data on chlorophyll is more limited and missing in smaller bodies of water, which is why a model based on data with this variable is based on fewer observations and fewer pseudoabsences. <br />
<br />
Information about which files are to be used as input in models for different species is gathered in the file "data.table.csv" which has the following format: <br />
<br />
<img src=images/table.input.jpg width=800 align=left> <br />
<br />
<br />
<br />
<br />
The script that loads data into the model can use different files for positive and negative findings is a inherited function from the previous project. In this project, positive and negative findings are read from the same file and the same pseudoabsences are used for all species, but it is possible, for example, to enter negative findings manually or have different "background models" for different groups of organisms by choosing different pseudoabsence files. <br />
<br />
The selection of points in the pseudoabsence model has not been filtered to remove duplicates with the same coordinates before sampling. This means that areas with many duplicates may be overrepresented. Filtering of duplicates was not done until the last run, otherwise it would of course have been done here as well. <br />
<br />
The "model" used to handle differences in sampling effort is one of the things that can be further developed in the next project. <br />
<br />

## Download and formatting of environmental data layers

We downloaded and formatted the following environmental data layers from https://neo.gsfc.nasa.gov/ and generated consistent continental data layers covering land and sea at 5 arcmin resolution for the following variables: Due to problems with downloading the larger files with data in °C and mg/m3 respectively and nonresponsive technical support we instead used the smaller image files with data as 1 byte/pixel in arbitrary units (AU). The classifier used, Random Forests, in a nonlinear model and sets cutoff based on the rank of values. Thus, scaling of data will not change the model and this time-saving shortcut does not impact the predictions. <br />
<br /> 
Sea surface temperature (SST)
* Annual mean SST in AU 
* Annual max SST in AU
* Annual min SST in AU
* Annual amplitude SST in AU <br />
<br />

Chlorofyll (Chla) 
* mean Chla of most productive month in AU 
* min Chla of most productive month in in AU
* max Chla of most productive month in in AU 
* amplitude Chla of most productive month in in AU <br />
<br />

Sea surface salinity (SSC)
* Annual mean SSS in PSU <br />
<br />

Layers for SST, Chla, and SSC are processed in the script Prepare environmental layers.HAV2022.r”. We used the source https://neo.gsfc.nasa.gov/archive/geotiff/ for data on Chlorophyll and SST respectively. Files with monthly data were downloaded from the following links:
https://neo.gsfc.nasa.gov/archive/geotiff/MY1DMM_CHLORA/  for chlorophyll 
https://neo.gsfc.nasa.gov/archive/geotiff/MYD28M/ for SST<br />
<br />
To create a salinity layer we took data for salinity in the oceans from Bio-Oracle (https://www.bio-oracle.org/). For simplicity, water not included in this data layer was assumed to be fresh water and the salinity was set to the limit for this at 0.05%. In practice by first setting all NA pixels to this value and then setting all pixels for which SST is missing to NA. <br />
<br />
 
## Model test and projection
 
Modellling was performed using the script masterscript.HAV2022.r with help functions library script SEanalytics.functionsNEW.r. The script includes the following steps <br />
<br />

* Filepaths, filenames, an suffixes
* Read present, absent and pseudoabsence points and extract environmental data
* Preparing iterations 
* Estimate weight of predictors with MCMC algortihm
* Make plots of variables’ weight and correlation with species observations
* Train random forest models
* Extract C50 rules 
* Calculate and plot ROC curves from cross validation 
* Spatial prediction of species presence and plotting individual maps 
* Plot map stacks with average probability
* Plot maps that combine cumulative average probability for all species with traffic layers <br />
<br />

### Filepaths, filenames, and suffixes
For the most important folders, three variants are created with different "suffixes" for results with and without the chlora variable. <br />
<br />

### Read present, absent and pseudoabsence points and extract environmental data
The filtered data from GBIF is loaded. Duplicates are then filtered out if they have the same coordinates and occurrence status. Next, environmental data is extracted from the raster stack that has the "correct" suffix. Observations without complete environmental data are filtered out. In the loop, a table is also created with statistics for each species on the number of positive and negative findings, including pseudoabsences. Species with fewer than 5 unique, complete, findings are filtered out. <br />
<br />

### Preparing iterations 
The analysis is performed as 5 replicates of a 5x5 cross-validation. For each replicate data are permuted. Then, positive and negative findings are sampled, individually, to belong to one of five possible (other values at the CV level are possible) sets. When an observation is assigned an iteration it means that it will be part of the test set during that iteration of the cross-validation. <br />
<br />

Since many points are close to each other and would lead to overestimation of predictive power if one allowed nearby points to be included in both the test and training sets, the coordinates are rounded as follows: <br />
<br />

* lonmin <-   10*floor(my.data$Lon/10)
* lonmax <-   10*ceiling(my.data$Lon/10)
* latmin <-   10*floor(2*my.data$Lat/10)/2
* latmax <-   10*ceiling(2*my.data$Lat/10)/2 <br />
<br />

and each observation is given an observation area whose name is given by the above. When you then sample data for training and test sets, observation areas are sampled, not individual observations. For some species, all findings ended up in a few observation areas. <br />
<br />

In order to carry out a 5x cross-validation, there must be findings in at least 5 different areas, otherwise the algorithm crashes. Since a meaningful estimate of predictive power, and also confidence in the maps, is doubtful if the positive findings come from fewer than five areas, such species are excluded. <br />
<br />

The size of the rasters and the very principle of grouping findings that may be considered too close to each other for meaningful ROC analysis can be debated. <br />
<br />

Information about the iterations is saved in a .rda file and is then used partly in the variable selection algorithm and partly when the Random Forest model is trained for cross-validation. <br />
<br />

### Estimate weight of predictors with MCMC algorithm
The weight of the variables is estimated with MCMC feature selection according to the same principle as in earlier model trials (Kruczyk et al 2012). Since the number of variables is small, however, we do not use the possibility of making a Random Forest model with only significant variables, but MCMC is only used to get a measure of the usefulness of the variables. The MCMC algorithm is run 25 times with different parts of the dataset to get a feel for how sensitive the weight of the variables (RI index) is to the selection of data. However, as the method is implemented, nearby points may be included in the same dataset and the weight of the variables may be affected by overtraining. This should be handled in the next project, for example by filtering out nearby finds. <br />
<br />

The MCMC algorithm is time-consuming and the 25 iterations are carried out on different nodes with the program package parallel{}. <br />
<br />

### Make plots of variables’ weight and correlation with species observations 
In part, a plot is made per species that shows the weight of the variables (RI index). A couple of random variables have been added to the model to prevent it from crashing. In the plot, a horizontal line has been inserted corresponding to the RI index that the model considers to be "better than chance". In addition, plots are made for each variable where x is the measured value of the predictor variable divided into 10 equally sized "bins" and the axis is the proportion of positive findings given that x lies in this bin. <br />
<br />
 
### Train random forest models
A model is trained for each sub-dataset according to the file created under “prepare iterations”. For these models, only the "probability that the observation has the status present" is saved for each observation. In addition, a model is trained with the entire input data that is saved and used for spatial prediction. <br />
<br />

In this implementation, iterations of Random Forests are run sequentially. The time gain of running in parallel on different nodes has not been that great, as the number of variables is small, and it has not been worth the time to modify the script for parallel training of models. However, this could be done easily in the future. <br />
<br />

### Extract C50 rules 
Test of a simple algorithm to create simple "rules" that describe the relationships between predictor variables and occurrence status. The output of this is saved as a text file for each species. The algorithm can be seen as an attempt to explain what happens inside the "decision trees" even if it is a different decision tree algorithm. The function has been tested without any optimization whatsoever, but it might be interesting to compare with the plots made above and as a "demo" of a track to follow up in future projects. <br />
<br />
 
### Calculate and plot ROC curves from cross validation 
ROC curves are calculated based on the Random Forest simulations above. The curves are based on the probability predicted for each observation when it constituted the test data, and thus with a model not based on points in this rectangle. Since the Random Forest algorithm was run in 5 repetitions of a 5x cross-validation, 5 ROC curves with different AUCs are obtained. All five curves are shown as a line in the plot and the average curve as a green polygon. How much difference there is between the AUC values from different repetitions gives an indication of how sensitive the model is to the data. <br />
<br />
 
### Spatial prediction of species presence and plotting individual maps (Sweden and global)
The Random Forest model trained with the entire dataset is used for prediction with a raster stack corresponding to the data used to train the model. The results are first saved as .rda files and, in a later step, as GeoTIFF. There is room here to make the script more uniform and skip the first step. <br />
<br />
 
### Plot map stacks with average probability
For each experiment, a raster stack is made with predicted probability from all species for which the model succeeded in creating such a map. An average value is calculated with the mean() function. There are also maps with logarithmic probability and weighted logarithmic probability. This is a legacy from the previous project by Bergkvist et al (2020). <br />
<br />
 
### Plot maps that combine cumulative average probability for all species with traffic layers
Raster layers with traffic data from 2016 that were used in the previous project (Bergkvist et al 2020) are read in and then the raster map corresponding to average probability is cropped and resampled to the same extent and resolution as the traffic map. The area being analyzed is given by the extent of the map used by Bergkvist et al (2020). <br />
<br />

A composite plot is made with basically the same code as in the previous project. The plot has four panels with <br />
<br />
* Average probability for all species given type of model (with or without chlorophyll and without chlorophyll but only points where chlorophyll data is available)
* Traffic data from 2016. The difference to the previous analysis is that the ceiling for traffic intensity was set to 1000 instead of 10000 (values above the ceiling are set to the ceiling). The reason for truncating high values is that otherwise you don't see the traffic lanes in the Gulf of Bothnia, which are much less intensive than the major waterways in the southern Baltic Sea.
* A map where the traffic data is superimposed on the model results with a different color.
* A map where traffic data and model results have been added. During the addition, a scale parameter is used which is set arbitrarily so that the traffic lane will appear.
 
### References
Kruczyk, M., H. Zetterberg, O. Hansson, S. Rolstad, L. Minthon, A. Wallin, K. Blennow, J. Komorowski and M. G. Andersson (2012). "Monte Carlo feature selection and rule-based models to predict Alzheimer’s disease in mild cognitive impairment." Journal of Neural Transmission 119(7): 821-831. <br />
<br />

Bergkvist J, Magnusson M, Obst M, Sundberg P, Andersson G (2020) Provtagningsdesign för övervakning av främmande arter. Övervakning i marin miljö. Havs- och vattenmyndighetens rapport 2020:22. ISBN 978-91-88727-86-2 <br />
<br />
