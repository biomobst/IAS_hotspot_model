library(dplyr)
library(readr)  
library(rgbif) # for occ_download

#1

# fill in your gbif.org credentials 
user <- "gunnarsva" # your gbif.org username 
pwd <- "Sak1arbiff" # your gbif.org password
email <- "gunnar.andersson@sva.se" # your email 

path <- "~/Dokument/Projekt/HAV2025/"
#file_url <- paste(path,"Data2022/valda.arter.5.csv",sep ="" )
file_url <- paste(path,"data/species.data/NIS_list_combined_Mar2025_v2.csv",sep ="" )


########################################
gbif_taxon_keys <- 
  readr::read_delim(file_url, delim =",",na = c("", "NA"), comment = "",   col_names = TRUE,skip_empty_rows = TRUE)%>%
  pull("Taxon name") %>% # use fewer names if you want to just test 
  name_backbone_checklist()  %>% # match to backbone
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) # get the gbif taxonkeys

gbif_taxon_keys <- gbif_taxon_keys

occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

###########################
d <- occ_download_get('0010903-240202131308920') %>%
  occ_download_import()