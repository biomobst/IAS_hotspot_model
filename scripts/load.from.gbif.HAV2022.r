library(dplyr)
library(readr)  
library(rgbif) # for occ_download

# fill in your gbif.org credentials 
user <- "username" # your gbif.org username 
pwd <- "password" # your gbif.org password
email <- "nn@host.country" # your email 


##################### Read list of species from csv file 
path <- "~/path/to/project/" # The path to the root of your project
file_url <- paste(path,"Data2022/valda.arter.4.csv",sep ="" ) # Filename for your list of species
#file_url <- paste(path,"Data2022/valda.arter.4.new.only.csv",sep ="" )# alternative filename


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

# A link to the result will be sent to your email when request is processed
###########################