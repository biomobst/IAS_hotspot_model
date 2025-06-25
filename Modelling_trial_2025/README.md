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


 
### References

Daraghmeh N, Exter K, Pagnier J, et al (2024) A long-term ecological research data set from the marine genetic monitoring programme ARMS-MBON 2018-2020. Molecular Ecology Resources. doi: https://doi.org/10.1111/1755-0998.14073

Pagnier J, Daraghmeh N, Obst M (2024) Using the long-term genetic monitoring network ARMS-MBON to detect marine non-indigenous species along the European coasts. Biological Invasions 27, 77 (2025). https://doi.org/10.1007/s10530-024-03503-2

Sundberg P, Obst M, Panova M (2024) DNA-baserad övervakning av arter i akvatisk miljö - verifiering och tillämpning. Naturvårdsverket rapport. pp 44. Rapport 7157. ISBN 978-91-620-7157-8
