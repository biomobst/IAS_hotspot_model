
Below is the explanation of the various functions and commands of the scripts used to perform Ecological niche modelling for invasive species (by Gunnar Andersson, National Veterinary Institute, SVA, Sweden, gunnar.andersson@sva.se)

<img src="/images/IAS_image4.png" width=550>

*Schematic representation of the file structure used by the script*

## The following files are included

**SEanalytics.functions.r**

The file contains all functions called in the other scripts.

**read.and.resample.raster.r**

The first script to run. 

The script extracts the zipfiles with the rasterlayers in a temporary directory. The compareRaster function is run to verify that all rasters have the same extent and projection, and thus can be included in a rasterstack. The rasters are cropped to match the region of interest and saved in temporary directory. Finally the rasters are aligned in a “rasterstack” which is saved as a .tif file. In the example stacks are prepared with a Global and a European limit.

The Global rasterstack will be used to extract environmental data for the traiingset whereas the limited European stack is used for predictions.

**extract.data.r**

Within script "read.and.resample.raster.r"

A test script to extract environmental data for one species from the stored stack.

**data.table.csv**

This file contains information for each species on which files to include for presence, absence and pseudoabsence data

**IASmodelling_speciesList_20200107.csv**

Contains information about classification of species. Information herein is used to prepare composite rasterlayers and identifying hotspots

**invasive.data.masterscript.r**

The first version of the script. A simpler implementation without cross validation

**invasive.data.masterscript.with CV**

This script performs all steps in the analysis, from extraction of data to preparation of plots. 

The analysis in performed in the following steps

**read.and.extract** 

Within script “SEanalytics.functions.r”

Read in present absent and pseudoabsent points, convert these points to spatial coordinates and extract environmental variables from rasterstack

Data points with errors are deleted and the data table, with environmental variables is saves as a .csv file

**split.data**

Within script “SEanalytics.functions.r”

The cross validation scheme is prepared. M repetitions of a n-fld cross validation. The cross validation is prepared in such a way that observations from the same ICES-statistical rectangle will end up in the same set. This is a precaution to avoid overestimating the predictive performance

**run.random.forests**

Within script “SEanalytics.functions.r”

The function “run random forests” is preparing the data and then the method is execuded in the sub-function RF-process. First the cross-validation experiments is preformed and predicted class probability for each observation is stored. (but not the modell). Finally a model is training with the full dataset and stored for future use. In this case class predictions for the training set are stored.
n.b. The code contains some “artifacts” that are left from an earlier implementation using parallel processing. The result is stored as two .rda files. One with cross validation  results the other with the full model

**ROC curve**

Within script “SEanalytics.functions.r”

Results are loaded from the cross validation results. The predicted class probability for each observation from the cross validation are sorted. FP TP FN and TN are calculated with different cutoffs. 

One average ROC curve is computed with results from all repetition of the cross validation as well as individual curves for each repetition. The meaning of repeating the cross validation is to check that the performance estimate is stable,  and does not vary a lot depending on how the data is split. 

**plot.maps**

Within script “SEanalytics.functions.r”

The rasterstack representing the area of interest is reloaded in memory. The function “plot maps” will run the random forests model on the raster. It returns the result as a raster with is then stored as an rda file. However, the function also makes some plots. The plots in the “plot maps” function could be inactivated to save time. However, it was useful to get the plots to see that the algorithm worked as intended.

For plotting the data different color plaettes are used using functions divpalette and seqpalette. Cutoffs are defined manually so that the color scale on the side will get a reasonable appearance and so that all species will be plotted with the same colorscale. In the end of the colorscale I add extra color for values outside the scale (that is land..)

Maps ar plotted with linear probabilities and log10 probabilities. To avoid na´s I add a small number to the probabilities before logging. 

After transformation the predicted rasters are saved as geotiff rasters.

**Combine layers**

Within script "invasive.data.masterscript"

A number of computations are preformed to combine results from several species. Information about what species belong to the same class (ecological or taxonomic) and their expected risk “samlat riskutfall” is read from a .csv file and stored in a table “combine.table”. The resulting combined layers are stored as geotiff.

**Plot selected layers**

Within script "invasive.data.masterscript.r"

The combined layers, stored as geotiff, are read in and plotted as png files.  Some computatons are done to obtain a reasonable colorscale so that plots are comparable and that the color-legend at the side will get a nice appearance (without breaks are strange numbers.)

Plots are made with either the entire predicted area (limited by the rasterstac) and for a smaller area defined by xlim and ylim.  The plots are overlayed with a map showing anministrative areas of Sweden and Europe. “shape 2”. 
