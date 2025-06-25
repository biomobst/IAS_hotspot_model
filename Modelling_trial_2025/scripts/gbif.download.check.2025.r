##require(rgdal)
require(CoordinateCleaner)
require(speciesgeocodeR)# installed via zip-tar frpm CRAN
require(sf)
require(raster)
#path_home
path <- "~/Dokument/Projekt/HAV2025/"
#file <- "Data2022/Corbicula fluminea GBIF/occurrence.txt"

#file <- "Data2022/GBIF alla arter 0429558-210914110416597/GBIF alla arter 0429558-210914110416597.csv"
#download 1, freshwater file <- "data/speciesdata/0010903-240202131308920.csv"
#file <- "data/speciesdata/0013730-240202131308920.csv"

#file <- "data/speciesdata/0011294-240216155721649.csv" #2 arter 22a feb
file <- "data/species.data/0018230-250310093411724.csv"# Arter från Matthias Mars 2024

input <- read.delim(paste(path,file,sep=""),header =T,sep="\t", 
                    na = c("", "NA"))#[ 1:1000,]
to.use <- complete.cases(input[,c("decimalLatitude" ,"decimalLongitude")])

input <- input[to.use,]

names(input)
#selcol <- c(1,3,12,16,62,68,70,71,71,73,76,77,84,108,108,109,115,116,122,137,138,139,175,189)
#head(input[,selcol])


################################# ################################# #################################
################################# plot  ################################# #################################
unique(input$species)
factor <- as.numeric(as.factor(input$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)

unique(fac1)
#names(input)[selcol]

#unique(input[,84])# occurance status
#which(input[,84]== "ABSENT")
#which(input[,84]== "PRESENT")

#unique(input[,112])# habitat

#unique(input[,115])# sampling effort
#unique(input[,137])# locationremarks


 #<- readOGR(dsn=paste(path,"data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp",sep=""),layer="CNTR_RG_01M_2020_4326")

world1 <- st_read(paste(path,"data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp",sep=""),layer="CNTR_RG_01M_2020_4326")
#plot(world1)
#If you only want the outlines, i.e. only the geometry, you need to plot the geometry part of the sf object. Try plot(st_geometry(NJ_Map_Road)) or plot(NJ_Map_Road$geometry), both should work.
plot(world1$geometry)

#xlim <- c(-180,180)
#ylim <- c(-60, 84)
par(mfrow=c(1,1))
xlim <- c(-100,20)
ylim <- c(8, 68)
#jpeg("testplot2.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
jpeg(paste(path,"testplots/testplot.mars 2025 smaller.jpg",sep=""), width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(input$"decimalLongitude",input$"decimalLatitude", col = fac2, pch = fac1)
dev.off()

################################# ################################# #################################
######################## clean coordinates #################################
################################# ################################# #################################
summary(input)
names(input)[c(22:23)]
cleanput <- input[-c(1),]
cleanput[,23] <- as.numeric(cleanput[,23])
names(cleanput)[c(22:23)] <- c("decimalLatitude" , "decimalLongitude")
hist(cleanput[,22])
hist(cleanput[,23])
which(!is.numeric(cleanput[,22]))
cleanput[,22]<- as.numeric(cleanput[,22])

head(cleanput)
cleanput <- clean_coordinates(x = cleanput)
summary(cleanput)

jpeg(paste(path,"testplots/plot.cleanput.result.mars 2025.jpg", sep=""), width = 1000,height=1000, pointsize = 10)
plot(cleanput)
dev.off()

save(cleanput, file = paste(path, "data/species.data/clean.coordinates.output.mars 2025.rda", sep=""))
# lines below used to merge datasets
#cleanput2 <- cleanput
#load(paste(path, "data/species.data/clean.coordinates.output.mars 2025.rda", sep=""))

#Merge two downloads
#cleanput <- rbind(cleanput,cleanput2)
#################################### plotta på världskarta igen ################################################
factor <- as.numeric(as.factor(cleanput$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)

unique(fac1)
jpeg(paste(path,"testplots/testplot.cleanput.mars2025.jpg",sep=""), width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$"decimalLongitude",cleanput$"decimalLatitude", col = fac2, pch = fac1)
dev.off()


####################################artvisa plottar ################################################
all.species <- unique(cleanput$species)
xlim <- c(-180,180)
ylim <- c(-40, 80)

#for(s in all.species[c(8,9)]){
  for(s in all.species){
    
jpeg(paste(path,"speciesplots/",s,"testplot.cleanput.mar2024.jpg",sep=""), width = 18*( xlim[2] - xlim[1]),height = 18*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$"decimalLongitude"[which(cleanput$species == s)],cleanput$"decimalLatitude"[which(cleanput$species == s)], 
       col = ifelse(cleanput$occurrenceStatus == "PRESENT",2,4),
       pch = ifelse(cleanput$occurrenceStatus == "PRESENT","+","o"), cex=4)
text(-170, -10, paste(s, "number of findings=", length(which(cleanput$species == s))), cex=20)
p <-length(which(cleanput$species == s & cleanput$occurrenceStatus == "PRESENT"))
a <- length(which(cleanput$species == s & cleanput$occurrenceStatus == "ABSENT"))
text(-170, -20, paste("n presence =", p, ": n absence =", a), cex=20)

title(main = s, cex=20)
dev.off()
}

############################## check if filter worked  
head(cleanput)
unique(cleanput$basisOfRecord)
for( u in unique(cleanput$basisOfRecord)){
  n <- length(which(cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
#}
#[1] "HUMAN_OBSERVATION"   "PRESERVED_SPECIMEN"  "FOSSIL_SPECIMEN"     "OCCURRENCE"          "MATERIAL_CITATION"   "MATERIAL_SAMPLE"    
#[7] "OBSERVATION"         "MACHINE_OBSERVATION" "LIVING_SPECIMEN"

exclude.Basis <- c("MACHINE_OBSERVATION","PRESERVED_SPECIMEN","FOSSIL_SPECIMEN", 
             "MATERIAL_CITATION"  , "MATERIAL_SAMPLE" ,"LIVING_SPECIMEN")
exclude.Basis.pattern <- paste(exclude.Basis, collapse = "|")
exlude.index <- grep (exclude.Basis.pattern, cleanput$basisOfRecord)

filtered.cleanput <- cleanput[-exlude.index,]
unique(filtered.cleanput$basisOfRecord)
for( u in unique(filtered.cleanput$basisOfRecord)){
  n <- length(which(filtered.cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
save(filtered.cleanput, file = paste(path, "data/species.data/filtered.clean.coordinates.output.mars 2025.rda", sep=""))

#save(filtered.cleanput, file = "filtered.clean.coordinates. mars 2024.output.rda")
###################### speciesgeododeR
#load( paste(path, "data/species.data/filtered.clean.coordinates.output.mars 2025.rda", sep=""))
#####load( "filtered.clean.coordinates. download feb 15.output.rda")


########## filter out data outside marine environment
stackpath <- paste(path, c("data/Biooracle.download/rasterstacks/baseline"), sep="/")
filename <- paste(stackpath,"/","filled_layers_new.tif", sep="")
mystack <- raster(filename)
mask <- mean(mystack)
plot(mask)

observations <- unique(filtered.cleanput[,c("gbifID","decimalLongitude", "decimalLatitude", "occurrenceStatus")])
names(observations) <- c("ID","Lon","Lat","occurrenceStatus")

#coordinatsystem. WGS 84 decimal-latitude. epsg:4326
observations$Lon <- as.numeric(observations$Lon)
observations$Lat <- as.numeric(observations$Lat)
coord<-as.data.frame(observations[,c("Lon","Lat")])
names(coord) <- c("Lon","Lat")

observations1 <- observations
#coord <- coord[c(1:5),]
#coordinates(coord) <- ~Lon+Lat

head(observations)
observations2 <- as.data.frame(observations[,-c(2,3)])
#points<-SpatialPointsDataFrame(coord,
#                              observations2, proj4string=CRS("+init=epsg:4326")) # expression decapriated
points<-SpatialPointsDataFrame(coord,
                               observations2, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xlim <- c(extent(points)[1],extent(points)[2]) + c(-0.01, 0.01)
ylim <- c(extent(points)[3],extent(points)[4]) + c(-0.01, 0.01)

extent(mask)

xlim2 <- extent(mask)[1:2]
ylim2 <- extent(mask)[3:4]

points2<-extract(mask, points, sp=TRUE)#t

head(points2)
#plot(mask)
plot(world1$geometry)
points(points, col=ifelse(is.na(points2$layer),2,3))

head(points2)
names(points2)


is.marine <- sapply(filtered.cleanput$gbifID, function(i)
  points2$layer[which( points2$ID == i)]
)
filtered.cleanput.marine <- filtered.cleanput[!is.na(is.marine),]
save(filtered.cleanput.marine, file = paste(path, "data/species.data/filtered.clean.marine.coordinates.output.mars 2025.rda", sep=""))
##############################################################
##################################################################
all.species <- unique(filtered.cleanput.marine$species)

#for(s in all.species){
#   s <- all.species[18]
#  temp <- filtered.cleanput[ which(filtered.cleanput$species == s & filtered.cleanput$occurrenceStatus == "PRESENT"),
#                    c("occurrenceID", "decimallongitude","decimallatitude")]
#  geocode1 <- SpeciesGeoCoder(temp, world1, areanames = "CNTR_NAME")

head(filtered.cleanput.marine)
#unique(filtered.cleanput$occurrenceStatus)
# s <- all.species[18]
speciespath <- paste(path, "/data/Indata2025/speciesIndata2025/", sep ="")
#prepare input file
selected.species <- read.csv2(file_url, sep=",")
cathegories <- unique(selected.species$category)
#my.cathegory <- cathegories[1]
tab <- c()
for(s in all.species){
  
  id <- which(!is.na(match(selected.species$Taxon.name,s)))
  if(length(id)<1){
    string <- strsplit(s," ")[[1]][1]
    id <- grep(string,selected.species$Taxon.name)
    comment <- selected.species$Taxon.name[id]
  }else{comment <- "names match"}
  cat <- selected.species$category[id]
  
  temp <- filtered.cleanput.marine[ which(filtered.cleanput$species == s ),
                             c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                               "depth", "depthAccuracy","eventDate")]
  #fix names of deciomal coordinates. capital L
  names(temp)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                  "depth", "depthAccuracy","eventDate")
  #head(temp)
  write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F)
  print(paste("wrote: ", speciespath,s,".csv", sep =""))
  print(paste("no positives: ",length(which(temp$occurrenceStatus == "PRESENT")), sep =""))
  print(paste("no negatives: ",length(which(temp$occurrenceStatus == "ABSENT")), sep =""))
  my.filename <- paste(s,".csv",sep="")
  my.pseudoname <- paste("pseudoabsences.marine.excludebox",cat,".csv",sep="")
  tab <- rbind(tab, c(s, my.filename,my.filename,my.pseudoname,
                      length(which(temp$occurrenceStatus == "PRESENT")),length(which(temp$occurrenceStatus == "ABSENT")),comment))
  


}
colnames(tab)<- c("species","present.data","absence.data","pseudoabsence.data","n.present","n.absent","comment on Taxon")
write.csv2(tab, file=paste(path, "/data/species.data/data.table.mars2025.csv",sep=""))
###################################################################
###################################################################
filtered.cleanput<-filtered.cleanput.marine
#define a box to exclude pseudoabsences
boxxlim=c(12,30)
boxylim = c(50,66)
xlim <- c(-180,180)
ylim <- c(-60, 84)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
lines( boxxlim[c(1,2,2,1,1)], boxylim[c(1,1,2,2,1)]) 

excludebox <- intersect(
  which(filtered.cleanput$decimalLongitude > boxxlim[1] & filtered.cleanput$decimalLongitude < boxxlim[2]),
  which(filtered.cleanput$decimalLatitude > boxylim[1] & filtered.cleanput$decimalLatitude < boxylim[2])
)
#plot excluded points

points(filtered.cleanput$decimalLongitude[excludebox], filtered.cleanput$decimalLatitude[excludebox], col=2)
#plot included points
points(filtered.cleanput$decimalLongitude[-excludebox], filtered.cleanput$decimalLatitude[-excludebox], col=3)

filtered.cleanput.unbox <- filtered.cleanput[-excludebox,]

##############################################
######### define groups of species      ######

file_url <- paste(path,"data/species.data/NIS_list_combined_Mar2025_v2.csv",sep ="" )
selected.species <- read.csv2(file_url, sep=",")
cathegories <- unique(selected.species$category)

#my.cathegory <- cathegories[1]
for (my.cathegory in cathegories[-5]){

my.species.list <- selected.species$Taxon.name[which(selected.species$category == my.cathegory)]

#get entries matching 
filtered.cleanput.subset <- filtered.cleanput.unbox[which(!is.na(match(filtered.cleanput.unbox$species,my.species.list))),]

locationsamples <- sample(1:length(filtered.cleanput.subset$gbifID), 1000, replace =F)

#locationsamples <- sample(1:length(filtered.cleanput.unbox$gbifID), 10000, replace =F)
pseudoabsences <- filtered.cleanput.subset[ locationsamples,
                                           c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                                             "depth", "depthAccuracy","eventDate")]
#pseudoabsences <- filtered.cleanput.unbox[ locationsamples,
#                                     c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
#                                       "depth", "depthAccuracy","eventDate")]
pseudoabsences$gbifID <- paste("pseudo",seq(1:1000), sep ="")
pseudoabsences$species <- NA
pseudoabsences$occurrenceStatus <- "ABSENT"
pseudoabsences$coordinateUncertaintyInMeters <- NA
pseudoabsences$depthAccuracy <- NA
pseudoabsences$eventDate <- NA
head(pseudoabsences)
names(pseudoabsences)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                "depth", "depthAccuracy","eventDate")

xlim <- c(-180,180)
ylim <- c(-60, 84)

#boxxlim=c(12,30)
#boxylim = c(50,66)
jpeg(paste(path,"/speciesplots/","testplot.pseudoabsences.",my.cathegory,".jpg",sep=""), width = 18*( xlim[2] - xlim[1]),height = 18*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
lines( boxxlim[c(1,2,2,1,1)], boxylim[c(1,1,2,2,1)]) 
points(pseudoabsences$"decimalLongitude",pseudoabsences$"decimalLatitude", col = "red", pch = "*", cex = 5)
dev.off()

write.csv(pseudoabsences,file = paste(speciespath,"pseudoabsences.marine.excludebox",my.cathegory,".csv", sep =""), row.names = F)
print(paste("wrote: ", speciespath,"/pseudoabsences.marine.excludebox",my.cathegory,".csv", sep =""))


}#end cathegory
##################xlim <- c(-100,20)
filtered.cleanput <- filtered.cleanput[ ,
                           c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                             "depth", "depthAccuracy","eventDate")]
#fix names of deciomal coordinates. capital L
#names(filtered.cleanput)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
#                "depth", "depthAccuracy","eventDate")

world1 <- st_read(paste(path,"data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp",sep=""),layer="CNTR_RG_01M_2020_4326")
#world1 <- readOGR(dsn=paste(path,"Data2022/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp",sep=""),layer="CNTR_RG_01M_2020_4326")

factor <- as.numeric(as.factor(filtered.cleanput$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)
xlim <- c(-110,40)
ylim <- c(-2, 68)#jpeg("testplot2.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
jpeg(paste(path,"speciesplots/","testplot.CLEANPUT.SMALLER.world.mar2025.jpg", sep=""), width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(filtered.cleanput$"decimalLongitude",filtered.cleanput$"decimalLatitude", col = fac2, pch = fac1)
dev.off()

#########################
