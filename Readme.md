{\rtf1\ansi\ansicpg1252\cocoartf2513
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Calibri;\f1\fswiss\fcharset0 ArialMT;\f2\fnil\fcharset0 Calibri-Bold;
\f3\ftech\fcharset77 Symbol;\f4\fnil\fcharset0 LucidaGrande;\f5\froman\fcharset0 TimesNewRomanPSMT;
\f6\fnil\fcharset0 Calibri-Italic;\f7\fmodern\fcharset0 CourierNewPSMT;\f8\fnil\fcharset0 Calibri-Light;
}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue255;\red24\green40\blue80;
}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c0\c0\c100000;\cssrgb\c12157\c21569\c38824;
}
\paperw11900\paperh16840\margl1440\margr1440\vieww14400\viewh8160\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs60 \cf2 \expnd0\expndtw0\kerning0
Documentation of hotspot model for the Gulf of Bothnia\
\pard\pardeftab720\partightenfactor0

\f1\fs45\fsmilli22667 \cf2 \
Summary \

\f0\fs32 \
This folder documents the a customisation of the original IAS hotspot model with the purpose to predict suitable habitat in the northern Baltic Sea for invasive species originating from freshwater systems, i.e. lakes and rivers. The model can be used to estimate the risk of limnic species spreading into or via the Gulf of Bothnia. We analyzed 34 species from different taxonomic groups and were able to produce individual models for 26 species, which you can find in this repo. The maps show that for many species found in northern European freshwater systems there are also suitable habitats in the Gulf of Bothnia. There is thus a risk for freshwater species to spread into the Gulf of Bothnia and other parts of the northern Baltic Sea.  The analysis of the models' predictions in relation to waterways that can spread the species between lakes and the sea and between different parts of the northern Baltic Sea justify monitoring port areas with large estuaries in the Gulf of Bothnia.
\f1\fs45\fsmilli22667 \
\
Download data from GBIF and pseudoabsences
\f0\fs32 \
\'a0\
Species data are downloaded from the input file \'93valda.arter.4.csv\'94 with the script \'94load.from.gbif.r\'94. Taxon names and Taxon keys are used to download occurrences (occ_download). A link to the data for download is sent to your email.\
\'a0\
Data filtration takes place in the scriptet \'94gbif.download.check.r\'94 , first by excluding erroneous coordinates with the function clean_coordinates() and thereafter by inspecting observations field \'94basisOfRecord\'94 and excluding the following classes: "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "MATERIAL_CITATION"\'a0 , "MATERIAL_SAMPLE" , "LIVING_SPECIMEN"). The following \'a0classes were included: "HUMAN_OBSERVATION", "OCCURRENCE", "OBSERVATION".\
\'a0\
For each species the cleaned data are saved as.csv file with the following columns: \
"gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "depth", "depthAccuracy","eventDate"\
\'a0\
For each species a csv file is generated with the following format \
\pard\pardeftab720\fi1738\partightenfactor0
\cf2 write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F)\
\pard\pardeftab720\partightenfactor0
\cf2 \'a0\
In addition to the files, the script creates plots of finds on a world map, partly for each species individually and partly a plot with all species in different colors to give an image of where in the world finds of invasive species have been made.\
\'a0\
The last step is to create pseudoabsences for presence-absence modeling. A plot of all downloaded observations (testplot.cleanput.jpg) shows that observations are unevenly distributed globally and that one can assume that there is a greater chance that someone found a species in an area where other species are reported. \
\'a0\
From the downloaded data files, cleaned data for 10,000 observations are sampled and used as pseudoabsences under the assumption that the distribution of all findings in "cleanput" can be considered a reasonable model for "sampling effort". Of the 10,000 findings, many will, in later stages, turn out to be at points where environmental data is missing, but in total about 1,500 pseudoabsences were useful, which fairly well balanced the positive findings for most species.\
\'a0\
In the tables "species stats" and "species stats no chlora" there is species-specific information on the number of finds in GBIF that meet the criteria above, the number of finds with unique coordinates and the number of these for which environmental information is available. Data on chlorophyll is more limited and missing in smaller bodies of water, which is why a model based on data with this variable is based on fewer observations and fewer pseudoabsences. \
\'a0\
Information about which files are to be used as input in models for different species is gathered in the file "data.table.csv" which has the following format:\
\'a0\

\itap1\trowd \taflags1 \trgaph108\trleft-108 \trwWidth11822\trftsWidth3 \trbrdrt\brdrnil \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalb \clshdrawnil \clwWidth2460\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx2160
\clvertalb \clshdrawnil \clwWidth3026\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx4320
\clvertalb \clshdrawnil \clwWidth3025\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx6480
\clvertalb \clshdrawnil \clwWidth2565\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f2\b\fs29\fsmilli14667 \cf2 species
\f0\b0\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f2\b\fs29\fsmilli14667 \cf2 present.data
\f0\b0\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f2\b\fs29\fsmilli14667 \cf2 absence.data
\f0\b0\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f2\b\fs29\fsmilli14667 \cf2 pseudoabsence.data
\f0\b0\fs32 \cell \row

\itap1\trowd \taflags1 \trgaph108\trleft-108 \trwWidth11822\trftsWidth3 \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalb \clshdrawnil \clwWidth2460\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx2160
\clvertalb \clshdrawnil \clwWidth3026\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx4320
\clvertalb \clshdrawnil \clwWidth3025\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx6480
\clvertalb \clshdrawnil \clwWidth2565\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Elodea nuttallii
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Elodea nuttallii.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Elodea nuttallii.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 pseudoabsences.csv
\fs32 \cell \row

\itap1\trowd \taflags1 \trgaph108\trleft-108 \trwWidth11822\trftsWidth3 \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalb \clshdrawnil \clwWidth2460\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx2160
\clvertalb \clshdrawnil \clwWidth3026\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx4320
\clvertalb \clshdrawnil \clwWidth3025\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx6480
\clvertalb \clshdrawnil \clwWidth2565\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Nymphoides peltata
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Nymphoides peltata.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Nymphoides peltata.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 pseudoabsences.csv
\fs32 \cell \row

\itap1\trowd \taflags1 \trgaph108\trleft-108 \trwWidth11822\trftsWidth3 \trbrdrl\brdrnil \trbrdrt\brdrnil \trbrdrr\brdrnil 
\clvertalb \clshdrawnil \clwWidth2460\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx2160
\clvertalb \clshdrawnil \clwWidth3026\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx4320
\clvertalb \clshdrawnil \clwWidth3025\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx6480
\clvertalb \clshdrawnil \clwWidth2565\clftsWidth3 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadl93 \clpadr93 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Pistia stratiotes
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Pistia stratiotes.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 Pistia stratiotes.csv
\fs32 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0

\fs29\fsmilli14667 \cf2 pseudoabsences.csv
\fs32 \cell \lastrow\row
\pard\pardeftab720\partightenfactor0
\cf2 \'a0\
The script that loads data into the model can use different files for positive and negative findings is a inherited function from the previous project. In this project, positive and negative findings are read from the same file and the same pseudoabsences are used for all species, but it is possible, for example, to enter negative findings manually or have different "background models" for different groups of organisms by choosing different pseudoabsence files.\
\'a0\
The selection of points in the pseudoabsence model has not been filtered to remove duplicates with the same coordinates before sampling. This means that areas with many duplicates may be overrepresented. Filtering of duplicates was not done until the last run, otherwise it would of course have been done here as well.\
\'a0\
The "model" used to handle differences in sampling effort is one of the things that can be further developed in the next project.\
\'a0\
\pard\pardeftab720\sa160\partightenfactor0

\f1\fs45\fsmilli22667 \cf2 Download and formatting of environmental data layers\
\pard\pardeftab720\partightenfactor0

\f0\fs32 \cf2 We downloaded and formatted the following environmental data layers from {\field{\*\fldinst{HYPERLINK "https://neo.gsfc.nasa.gov/"}}{\fldrslt \cf3 \ul \ulc3 https://neo.gsfc.nasa.gov/}}\cf3 \ul \ulc3  \cf2 \ulnone and generated consistent continental data layers covering land and sea at 5 arcmin resolution for the following variables: Due to problems with downloading the larger files with data in 
\f3 \'a1
\f0 C and mg/m3 respectively and nonresponsive technical support we instead used the smaller image files with data as 1 byte/pixel in arbitrary units (AU). The classifier used, Random Forests, in a nonlinear model and sets cutoff based on the rank of values. Thus, scaling of data will not change the model and this time-saving shortcut does not impact the predictions. \
\pard\pardeftab720\li480\partightenfactor0
\cf3 \ul \'a0\cf2 \ulnone \
\pard\pardeftab720\li480\fi-480\partightenfactor0

\f4 \cf2 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f6\i\fs32 Sea surface temperature
\f0\i0  (SST)\
\pard\pardeftab720\li1440\fi-480\partightenfactor0

\f7 \cf2 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 Annual mean SST in AU \

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 Annual max SST in AU\

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 Annual min SST in AU\

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 Annual amplitude SST in AU\
\pard\pardeftab720\li480\fi-480\partightenfactor0

\f4 \cf2 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f6\i\fs32 Chlorofyll
\f0\i0  (Chla) \
\pard\pardeftab720\li1440\fi-480\partightenfactor0

\f7 \cf2 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 mean Chla of most productive month in AU \

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 min Chla of most productive month in in AU\

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 max Chla of most productive month in in AU \

\f7 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 amplitude Chla of most productive month in in AU\
\pard\pardeftab720\li480\fi-480\partightenfactor0

\f4 \cf2 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f6\i\fs32 Sea surface salinity (SSC)
\f0\i0 \
\pard\pardeftab720\li1440\fi-480\partightenfactor0

\f7 \cf2 o
\f5\fs18\fsmilli9333 \'a0\'a0 
\f0\fs32 Annual mean SSS in PSU\
\pard\pardeftab720\partightenfactor0
\cf2 \'a0\
Layers for SST, Chla, and SSC are processed in the script Prepare environmental layers.HAV2022.r\'94. We used the source {\field{\*\fldinst{HYPERLINK "https://neo.gsfc.nasa.gov/archive/geotiff/"}}{\fldrslt \cf3 \ul https://neo.gsfc.nasa.gov/archive/geotiff/}} for data on Chlorophyll and SST respectively. Files with monthly data were downloaded from the following links:\
\pard\pardeftab720\partightenfactor0
{\field{\*\fldinst{HYPERLINK "https://neo.gsfc.nasa.gov/archive/geotiff/MY1DMM_CHLORA/"}}{\fldrslt \cf3 \ul https://neo.gsfc.nasa.gov/archive/geotiff/MY1DMM_CHLORA/}}\'a0 for chlorophyll \
\pard\pardeftab720\partightenfactor0
{\field{\*\fldinst{HYPERLINK "https://neo.gsfc.nasa.gov/archive/geotiff/MYD28M/"}}{\fldrslt \cf3 \ul https://neo.gsfc.nasa.gov/archive/geotiff/MYD28M/}} \'a0for SST\
\'a0\
To create a salinity layer we took data for salinity in the oceans from Bio-Oracle ({\field{\*\fldinst{HYPERLINK "https://www.bio-oracle.org/"}}{\fldrslt \cf3 \ul https://www.bio-oracle.org/}}). For simplicity, water not included in this data layer was assumed to be fresh water and the salinity was set to the limit for this at 0.05%. In practice by first setting all NA pixels to this value and then setting all pixels for which SST is missing to NA.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f1\fs45\fsmilli22667 \cf2 Model test and projection
\f0\fs32 \
\'a0\
Modellling was performed using the script masterscript.HAV2022.r with help functions library script SEanalytics.functionsNEW.r\
\'a0\
The script includes the following steps\
\pard\pardeftab720\li480\fi-480\partightenfactor0

\f4 \cf2 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Filepaths, filenames, an suffixes\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Read present, absent and pseudoabsence points and extract environmental data\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Preparing iterations \

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Estimate weight of predictors with MCMC algortihm\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Make plots of variables\'92 weight and correlation with species observations\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Train random forest models\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Extract C50 rules \

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Calculate and plot ROC curves from cross validation \

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Spatial prediction of species presence and plotting individual maps \

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Plot map stacks with average probability\

\f4 \'b7
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Plot maps that combine cumulative average probability for all species with traffic layers\
\pard\pardeftab720\partightenfactor0
\cf2 \'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Filepaths, filenames, and suffixes
\f0 \cf2 \
For the most important folders, three variants are created with different "suffixes" for results with and without the chlora variable.\
\'a0\

\f8 \cf4 Read present, absent and pseudoabsence points and extract environmental data\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 The filtered data from GBIF is loaded. Duplicates are then filtered out if they have the same coordinates and occurrence status. Next, environmental data is extracted from the raster stack that has the "correct" suffix. Observations without complete environmental data are filtered out. In the loop, a table is also created with statistics for each species on the number of positive and negative findings, including pseudoabsences. Species with fewer than 5 unique, complete, findings are filtered out.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Preparing iterations \
\pard\pardeftab720\partightenfactor0

\f0 \cf2 The analysis is performed as 5 replicates of a 5x5 cross-validation. For each replicate data are permuted. Then, positive and negative findings are sampled, individually, to belong to one of five possible (other values at the CV level are possible) sets. When an observation is assigned an iteration it means that it will be part of the test set during that iteration of the cross-validation.\
\'a0\
Since many points are close to each other and would lead to overestimation of predictive power if one allowed nearby points to be included in both the test and training sets, the coordinates are rounded as follows:\
\'a0 lonmin <-\'a0\'a0 10*floor(my.data$Lon/10)\
\'a0 lonmax <-\'a0\'a0 10*ceiling(my.data$Lon/10)\
\'a0 latmin <-\'a0\'a0 10*floor(2*my.data$Lat/10)/2\
\'a0 latmax <-\'a0\'a0 10*ceiling(2*my.data$Lat/10)/2\
and each observation is given an observation area whose name is given by the above. When you then sample data for training and test sets, observation areas are sampled, not individual observations. For some species, all findings ended up in a few observation areas.\
\'a0\
In order to carry out a 5x cross-validation, there must be findings in at least 5 different areas, otherwise the algorithm crashes. Since a meaningful estimate of predictive power, and also confidence in the maps, is doubtful if the positive findings come from fewer than five areas, such species are excluded.\
\'a0\
The size of the rasters and the very principle of grouping findings that may be considered too close to each other for meaningful ROC analysis can be debated.\
\'a0\
Information about the iterations is saved in a .rda file and is then used partly in the variable selection algorithm and partly when the Random Forest model is trained for cross-validation.\
\pard\pardeftab720\partightenfactor0

\f2\b \cf2 \'a0
\f0\b0 \
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Estimate weight of predictors with MCMC algorithm\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 The weight of the variables is estimated with MCMC feature selection according to the same principle as in earlier model trials (Kruczyk et al 2012). Since the number of variables is small, however, we do not use the possibility of making a Random Forest model with only significant variables, but MCMC is only used to get a measure of the usefulness of the variables. The MCMC algorithm is run 25 times with different parts of the dataset to get a feel for how sensitive the weight of the variables (RI index) is to the selection of data. However, as the method is implemented, nearby points may be included in the same dataset and the weight of the variables may be affected by overtraining. This should be handled in the next project, for example by filtering out nearby finds.\
\'a0\
The MCMC algorithm is time-consuming and the 25 iterations are carried out on different nodes with the program package parallel\{\}.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Make plots of variables\'92 weight and correlation with species observations 
\f0 \cf2 \
In part, a plot is made per species that shows the weight of the variables (RI index). A couple of random variables have been added to the model to prevent it from crashing. In the plot, a horizontal line has been inserted corresponding to the RI index that the model considers to be "better than chance". In addition, plots are made for each variable where x is the measured value of the predictor variable divided into 10 equally sized "bins" and the axis is the proportion of positive findings given that x lies in this bin.\
\'a0\

\f8 \cf4 Train random forest models
\f0 \cf2 \
A model is trained for each sub-dataset according to the file created under \'93prepare iterations\'94. For these models, only the "probability that the observation has the status present" is saved for each observation. In addition, a model is trained with the entire input data that is saved and used for spatial prediction.\
\'a0\
In this implementation, iterations of Random Forests are run sequentially. The time gain of running in parallel on different nodes has not been that great, as the number of variables is small, and it has not been worth the time to modify the script for parallel training of models. However, this could be done easily in the future. \
\'a0\

\f8 \cf4 Extract C50 rules \
\pard\pardeftab720\partightenfactor0

\f0 \cf2 Test of a simple algorithm to create simple "rules" that describe the relationships between predictor variables and occurrence status. The output of this is saved as a text file for each species. The algorithm can be seen as an attempt to explain what happens inside the "decision trees" even if it is a different decision tree algorithm. The function has been tested without any optimization whatsoever, but it might be interesting to compare with the plots made above and as a "demo" of a track to follow up in future projects.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Calculate and plot ROC curves from cross validation \
\pard\pardeftab720\partightenfactor0

\f0 \cf2 ROC curves are calculated based on the Random Forest simulations above. The curves are based on the probability predicted for each observation when it constituted the test data, and thus with a model not based on points in this rectangle. Since the Random Forest algorithm was run in 5 repetitions of a 5x cross-validation, 5 ROC curves with different AUCs are obtained. All five curves are shown as a line in the plot and the average curve as a green polygon. How much difference there is between the AUC values from different repetitions gives an indication of how sensitive the model is to the data.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Spatial prediction of species presence and plotting individual maps (Sweden and global)\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 The Random Forest model trained with the entire dataset is used for prediction with a raster stack corresponding to the data used to train the model. The results are first saved as .rda files and, in a later step, as GeoTIFF. There is room here to make the script more uniform and skip the first step.\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Plot map stacks with average probability\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 For each experiment, a raster stack is made with predicted probability from all species for which the model succeeded in creating such a map. An average value is calculated with the mean() function. There are also maps with logarithmic probability and weighted logarithmic probability. This is a legacy from the previous project by Bergkvist et al (2020). \
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 Plot maps that combine cumulative average probability for all species with traffic layers\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 Raster layers with traffic data from 2016 that were used in the previous project (Bergkvist et al 2020) are read in and then the raster map corresponding to average probability is cropped and resampled to the same extent and resolution as the traffic map. The area being analyzed is given by the extent of the map used by Bergkvist et al (2020). \
\'a0\
A composite plot is made with basically the same code as in the previous project. The plot has four panels with\
\'a0\
\pard\pardeftab720\li960\fi-960\sl384\partightenfactor0

\f8 \cf2 \'95
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Average probability for all species given type of model (with or without chlorophyll and without chlorophyll but only points where chlorophyll data is available)\

\f8 \'95
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0 
\f0\fs32 Traffic data from 2016. The difference to the previous analysis is that the ceiling for traffic intensity was set to 1000 instead of 10000 (values above the ceiling are set to the ceiling). The reason for truncating high values is that otherwise you don't see the traffic lanes in the Gulf of Bothnia, which are much less intensive than the major waterways in the southern Baltic Sea.\

\f8 \'95
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0 
\f0\fs32 A map where the traffic data is superimposed on the model results with a different color.\

\f8 \'95
\f5\fs18\fsmilli9333 \'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0\'a0 
\f0\fs32 A map where traffic data and model results have been added. During the addition, a scale parameter is used which is set arbitrarily so that the traffic lane will appear.\
\pard\pardeftab720\partightenfactor0
\cf2 \'a0\
\'a0\
\pard\pardeftab720\partightenfactor0

\f8 \cf4 References\
\pard\pardeftab720\partightenfactor0

\f0 \cf2 Kruczyk, M., H. Zetterberg, O. Hansson, S. Rolstad, L. Minthon, A. Wallin, K. Blennow, J. Komorowski and M. G. Andersson (2012). "Monte Carlo feature selection and rule-based models to predict Alzheimer\'92s disease in mild cognitive impairment." Journal of Neural Transmission 119(7): 821-831.\
\'a0\
Bergkvist J, Magnusson M, Obst M, Sundberg P, Andersson G (2020) Provtagningsdesign f\'f6r \'f6vervakning av fr\'e4mmande arter. \'d6vervakning i marin milj\'f6. Havs- och vattenmyndighetens rapport 2020:22. ISBN 978-91-88727-86-2\
\'a0\
}