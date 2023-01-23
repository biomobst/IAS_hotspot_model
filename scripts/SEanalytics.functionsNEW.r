#functions
require(randomForest)

require(raster)
require(rgdal)
require(fBasics)
#############################################################
################### read and extract rasterdata  ############
#############################################################
read.and.extract <- function(data.table,species,stack, speciespath,stackpath,plotpath,outpath){
  #print(species)
  presence <- data.table$present.data[which(data.table$species == species)]
  #print(head(presence))
  absence <- data.table$present.data[which(data.table$species == species)]
  # print(head(absence))
  pseudoabsence <- data.table$pseudoabsence.data[which(data.table$species == species)]
  #pseudoabsence <- strsplit(pseudoabsence, ":")[[1]]
  #print(head(pseudoabsence))
  #
  present.data <-read.csv(paste(speciespath,presence,sep="/"), stringsAsFactors = F)
  present.data <- present.data[which(present.data$occurrenceStatus == "PRESENT"),]
  names(present.data)[1] <-"ID"
  present.data$occurrenceStatus <- "present"
  
  if(is.na(absence)){absence <- ""}
  if(!absence==""){
    absent.data <-read.csv(paste(speciespath,absence,sep="/"), stringsAsFactors = F)
    absent.data <- absent.data[which(absent.data$occurrenceStatus == "ABSENT"),]
    names(absent.data)[1] <-"ID"
    if(length(absent.data$ID) >0){  absent.data$occurrenceStatus <- "absent" }
    
    if(length(absent.data$occurrenceID >0)){
      names(absent.data)[1] <-"ID"
      print(paste("absences found",species))
    }else{print(paste("No absences in absence file",species))}
    
  }else{
    absent.data = present.data[FALSE,]
    
    print(paste("No absences",species))
  }
  
  pseudoabsence.data <- c()
  for(i in pseudoabsence){
    temp <- read.csv(paste(speciespath,i,sep="/"),
                     sep=",", stringsAsFactors = F)[,c("gbifID","decimalLongitude",
                                                       "decimalLatitude","occurrenceStatus")]
    
    names(temp) <- c("ID", "decimalLongitude","decimalLatitude","occurrenceStatus")
    if(length(pseudoabsence.data) >0){
      pseudoabsence.data <- rbind(pseudoabsence.data,temp) 
    }else{pseudoabsence.data <- temp}
  }
  pseudoabsence.data$occurrenceStatus <- "absent"
  head(present.data)
  #head(absent.data)
  
  points.pres <- present.data[,c("ID","decimalLongitude","decimalLatitude","occurrenceStatus")]
  points.abs <- absent.data[,c("ID","decimalLongitude","decimalLatitude","occurrenceStatus")]
  
  
  points.pseudo <- pseudoabsence.data[,c("ID","decimalLongitude","decimalLatitude","occurrenceStatus")]
  points.all <- rbind(points.pres,points.abs,points.pseudo)
  dim(unique(points.all[,c("decimalLongitude", "decimalLatitude", "occurrenceStatus")]) )
  dim(points.all[,c("decimalLongitude", "decimalLatitude", "occurrenceStatus")]) 
  names(points.all) <- c("ID","Lon","Lat","occurrenceStatus")
  
  npos <- length(which(points.all$occurrenceStatus == "present"))
  nabs <- length(which(points.all$occurrenceStatus == "absent"))
  
  duplicates <- which(duplicated(points.all[, c("Lon","Lat","occurrenceStatus")]))
  points.all<- points.all[-duplicates,]
  
  npos.unique <- length(which(points.all$occurrenceStatus == "present"))
  nabs.unique <- length(which(points.all$occurrenceStatus == "absent"))
  
  points.all$Lon <- as.numeric(points.all$Lon)
  points.all$Lat <- as.numeric(points.all$Lat)
  coord<-points.all[,c("Lon","Lat")]
  names(coord) <- c("Lon","Lat")
  remove <- which(!complete.cases(coord))
  if(length(remove >0)){
    coord <- coord[-remove,]
    points.all <- points.all[-remove,]
  }
  points<-SpatialPointsDataFrame(coord,
                                 points.all, proj4string=CRS("+init=epsg:4326"))
  
  
  
  
  ## load rasterdata
  # "globalStack.rda" or "europeStack.rda"
  #load(paste(stackpath,stack, sep="/"))
  # testplot
  # xlim <- c(extent(points)[1],extent(points)[2]) + c(-0.01, 0.01)
  # ylim <- c(extent(points)[3],extent(points)[4]) + c(-0.01, 0.01)
  # extent(var2)
  # xlim2 <- extent(var2)[1:2]
  # ylim2 <- extent(var2)[3:4]
  #plot(mean(var2, na.rm=T), xlim=xlim2, ylim=ylim2)
  # proj4string(points)
  # proj4string(var2)
  # extract raster data for all points
  #if(stack ==  "globalStack.rda"){
  #  points2<-extract(var, points, sp=TRUE)#t
  #  rm(var)
  #}
  #if(stack ==  "europeStack.rda"){
  #  points2<-extract(var2, points, sp=TRUE)#t
  #  rm(var2)
  #}
  points2<-extract(stack, points, sp=TRUE)
  head(points2)
  
  names(points2)
  # filter out incomplete points (e.g points at land or outside )
  points2 <- as.data.frame(points2)
  complete.points <- points2[which(complete.cases(points2)),]
  head(complete.points)
  names(complete.points)
  # repove the extra colums for "Lat" and "Lon" that was inserted when creating the dataframe
  colnames(complete.points)
  complete.points <- complete.points[,-match(c("Lon.1", "Lat.1"), colnames(complete.points)) ]
  head(complete.points)
  # rm(var2)
  npos.unique.complete <- length(which(complete.points$occurrenceStatus == "present"))
  nabs.unique.complete <- length(which(complete.points$occurrenceStatus == "absent"))
  
  return(list("complete.points" = complete.points , 
              "stats" = c(npos,nabs, npos.unique, nabs.unique, npos.unique.complete, nabs.unique.complete)))
}
#############################################################
################### Run random forests  ############

#############################################################
################### Run random forests  ############
##################################################################
run.random.forests <- function(species, selvar, indata.path,iterations.path){
  lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))
  # load the information about iterations and cross validation for species in question
  load(paste(iterations.path,"/occurance.iters_",species,".rda",sep=""))
  
  ## remove corrupt data
  my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)
  print(paste(species,paste(unique(my.data$occurrenceStatus) ))  )
  my.data$occurrenceStatus <- as.character(my.data$occurrenceStatus)
  present.synonyms <-  which(my.data$occurrenceStatus == "present"|my.data$occurrenceStatus == "Present"| my.data$occurrenceStatus == "established"| my.data$occurrenceStatus == "Established")
 my.data$occurrenceStatus[present.synonyms] <- "present"
 absent.synonyms <-  which(my.data$occurrenceStatus == "Absent"| my.data$occurrenceStatus == "")
 
 my.data$occurrenceStatus[absent.synonyms] <- "absent"
 remove <- which(!my.data$occurrenceStatus == "absent"  & !my.data$occurrenceStatus == "present")
 print(paste("remove",remove))
 if(length(remove > 0)){
   my.data <- my.data[-remove,]                
 }
 
 ###
  print(paste(species,paste(unique(my.data$occurrenceStatus) ))  )
  my.data$occurrenceStatus <- as.character(my.data$occurrenceStatus)
  #make process plan based on iterations file...
  
  nrep <- length(all.occurance.iters)
  CV.level <- length(all.occurance.iters[[1]])
  process <- seq(1,nrep*CV.level)
  repeats <- sort(rep(seq(1,nrep),CV.level))
  iters <- rep(seq(1,CV.level),nrep)
  #process.plan <- cbind(process,rep,iter)
  
  ##################################
  ### run random forests in cross validation
  #####################################
  RF.results <- list()
  for( i in 1:nrep){
    RF.results[[i]]<- list()
  }
  for(i in 1:(CV.level*nrep)){
  # i <- 1
    iter <- iters[i]
    rep <- repeats[i]
      
  process.ID <-NA
  my.control <- list()
 
  #we may want to include a feature selection step. this this option
   if(selvar == "all"){
     my.control[["sel.var"]] <- colnames(my.data)[-c(1:4)]  }
 print(paste(species, "rep",rep,"iter", iter))

  # pos <- sample(my.data$ID[which(my.data$occurrenceStatus == "present")],
    #            floor(length(my.data$ID[which(my.data$occurrenceStatus == "present")])*0.95) , replace = FALSE)
  #neg <- sample(my.data$ID[which(my.data$occurrenceStatus == "absent")],
   #             floor(length(my.data$ID[which(my.data$occurrenceStatus == "absent")])*0.95) , replace = FALSE)
 #  my.control[["train"]] <- c(pos,neg)

 ## get the correct train data based on the MCMC version
 iter.ID <- all.occurance.iters[[rep]][[iter]]
 #iter.nr <-  unlist(sapply(1:length(iter.ID),function(i)
 #  grep(paste("A",iter.ID[i],"A"),paste("A",.mydata$ID,"A"))
 #))    
 #which(my.data$ID %in% iter.ID)
 
.temptrain <- my.data$ID[- which(my.data$ID %in% iter.ID) ]
 ##
 my.control[["train"]] <-  unlist( .temptrain ) ###all.occurance.iters[[rep]][[iter]]
 
 # for(r in 1:5){print(length( all.occurance.iters[[rep]][[r]] ))}
  RF.results[[rep]][[iter]] <-  RF.process(.mydata =my.data,
                           .my.control =my.control,
                          .iter = iter, 
                          .process.ID = process.ID,
                          .store.model=FALSE)
  }
  ################################
  ### run RF with all data in trainingset
  ##################################
  iter <- 1
  process.ID <-NA
  my.control <- list()
  
  #we may want to include a feature selection step. this this option
  if(selvar == "all"){
    my.control[["sel.var"]] <- colnames(my.data)[-c(1:4)]  }
  print(paste(species, "iter", iter))
 # pos <- sample(my.data$ID[which(my.data$occurrenceStatus == "present")],
 #               floor(length(my.data$ID[which(my.data$occurrenceStatus == "present")])*0.95) , replace = FALSE)
  #neg <- sample(my.data$ID[which(my.data$occurrenceStatus == "absent")],
    #            floor(length(my.data$ID[which(my.data$occurrenceStatus == "absent")])*0.95) , replace = FALSE)
  
  my.control[["train"]] <- my.data$ID
  
  RF.results.alldata <- RF.process(.mydata =my.data,
                           .my.control <- my.control,
                           .iter <- iter,
                           .process.ID <- process.ID,
                           .store.model=TRUE)
  
  
  return(list("RF.results.CV" = RF.results, "RF.results.alldata" = RF.results.alldata))
}

#############################################################
################### Run random forests process  ############
##################################################################
### this is the actual process that executes randomRorests. It is called by another funtion that prepare data
# on classweights https://stackoverflow.com/questions/20251839/how-to-use-classwt-in-randomforest-of-r
# wn = sum(y="N")/length(y)
#wy = 1
#Then set classwt = c("N"=wn,"Y"=wy)
#Alternatively, you may want to use ranger package. This package offers flexible builds of random forests, and specifying class / sample weight is easy. ranger is also supported by caret package

RF.process = function(.mydata, .my.control, .iter, .process.ID, .store.model=FALSE){
  .sel.var <- .my.control[["sel.var"]]# to doublecheck...
  
  .train <- .my.control[["train"]]# to doublecheck...
  .train.index <- which(.mydata$ID %in% .train)
  #length(which(.mydata$ID %in% .train))
  #length(.mydata$ID)
  #length(.train.index)
  
  #print(head(my.data))
  #print(sel.var)
 # print(.train.index)
  
  .remove.col <- unlist( sapply(c("ID","Lon" ,"Lat" ,"occurrenceStatus"  ),function(i)
    grep(i,names(.mydata))
  ))
  .value <- as.factor(.mydata$occurrenceStatus)
  .d <- .mydata[,-.remove.col]
  .sel.var.index <- match(.sel.var,colnames(.d))
  
  #sum( .d[1,] )
  
  print("RandomForests selected variables" )
#  RF.selected <- try(randomForest(.d[.train.index,.sel.var.index],.value[.train.index], 
  #                 xtest = .d[-.train,.sel.var.index], ytest =  .value[-.train.index] ,             
   #                               prox=TRUE) )
  
  w.pres = sum(.value[.train.index]=="present")/length(.value[.train.index])
  w.abs = 1
  print(paste("w.pres=",w.pres,"wabs=",w.abs))
  RF.selected <- try(randomForest(.d[.train.index,.sel.var.index],.value[.train.index],
                                  prox=FALSE,classwt = c("absent"=w.abs,"present"=w.pres)) )
  
  success <- class(RF.selected ) != "try-error"
  
  
  print("Predict.selected")
 # temp.predict.sel <- c()
#  cutoffs <- c(seq(1:99)/100)
  
  #if trainingset is not entire dataset, predict testset, else predict trainingset.
  if( length(.mydata$ID)>length(.train) ){
    
 # for(j in 1:length(cutoffs)){
 #   cutoff1 <- cutoffs[j]
#
 #   temp.predict.sel<- rbind(temp.predict.sel,
 #                            predict(RF.selected,.d[-.train,.sel.var.index],cutoff=c(cutoff1, 1-cutoff1)))
 # }
    print("predict alt 1")
    
  probs.sel <- predict(RF.selected,.d[-.train.index,.sel.var.index],type="prob")
  print("predicted alt 1")
  
  }else{
    
   # for(j in 1:length(cutoffs)){
   #   cutoff1 <- cutoffs[j]
   #   
    #  temp.predict.sel<- rbind(temp.predict.sel,
   #                            predict(RF.selected,.d[,.sel.var.index],cutoff=c(cutoff1, 1-cutoff1)))
  #  }
    print("predict alt 2")
  probs.sel <- predict(RF.selected,.d[,.sel.var.index],type="prob")   
  print("predicted alt 2")
  
  }
  print("store.model.check")
  #response.sel <- predict(RF.selected,.d[-.train,.sel.var.index],type="response")
  if(.store.model==FALSE){
    RF.selected <- NA
    print("store.model.check F")}
  my.return <- list("RF.selected" = RF.selected, 
                     "predict.selected" = NA,
                     "probs.sel" = probs.sel)
  print("store.model")
  
  return(my.return)
}

#############################################################
################### predict and plot the maps  ############
##################################################################
### 
plot.maps <- function(species,indata.path,modelpath,plotpath, colors,brk){
  lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))
  my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)
  
  load(paste(modelpath,"/RF.model.and.predictions.eur.wt.",species,".rda",sep=""))
  model <- rf.output$RF.selected
  map<-predict(bar, model, type="prob")
  map <- 1-map # to get probability of presence
  
  map2 <- map
  map2[is.na(map2)]<- 1.005
  colors <- c(colors,"lightgrey")
  
  xlim=c(0,30)# for swe
  ylim = c(50,70)#for swe
  
  # linear Europe
  png(paste(plotpath,"/",species, ".lin.trainpoints.world.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(map2,col=colors, breaks = brk)
  plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  # points(my.data$Lon, my.data$Lat, col = ifelse(my.data$occurrenceStatus == "absent", 1,2),
  #        pch=ifelse(my.data$occurrenceStatus == "absent", ".","."))
  points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4,  pch=1, cex=0.05, lwd=0.05)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2,  pch=1, cex=0.05, lwd=0.05)
  
  dev.off()
  #### linear Sweden
  png(paste(plotpath,"/",species, ".lin.trainpoints.Swe.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(map2,col=colors, breaks = brk, xlim=xlim, ylim=ylim)
  plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  # points(my.data$Lon, my.data$Lat, col = ifelse(my.data$occurrenceStatus == "absent", 1,2),
  #        pch=ifelse(my.data$occurrenceStatus == "absent", ".","."))
  points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4,  pch=2, cex=0.3, lwd=0.25)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2,  pch=2, cex=0.3, lwd=0.25)
  
  dev.off()
  
 # png(paste(plotpath,"/",species, "B.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
 # plot(map2,col=colors, breaks = brk)
 # plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
 #  dev.off()
  #png(paste(plotpath,"/C.",species, ".png",sep=""),  width = 180, height = 180, units = "mm", res=6000)
  # plot(map,col=colors)
  # plot(shape.ecoregions,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
   
   #dev.off()
   
   return(map)
}
#############################################################
################### split data for cross validation  ############
##################################################################
### 
split.data <- function(species,indata.path,iterations.path){
  lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))
  my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)
  
  length(my.data$ID)
  length(unique(my.data$ID))
  
  #  names(my.data)
  #my.data$Lat
  #my.data$Lon
  #Define in what ICES statistical rectangle the sample is in
  #recatangles are 1 degree longitude, 0.5 degree latitude#
  lonmin <-   10*floor(my.data$Lon/10)
  lonmax <-   10*ceiling(my.data$Lon/10)
  latmin <-   10*floor(2*my.data$Lat/10)/2
  latmax <-   10*ceiling(2*my.data$Lat/10)/2
 # define in which ICES statistical rectangle the datapoint is
 # lonmin <-   1*floor(my.data$Lon/1)
#  lonmax <-   1*ceiling(my.data$Lon/1)
 # latmin <-   1*floor(2*my.data$Lat/1)/2
#  latmax <-   1*ceiling(2*my.data$Lat/1)/2
  #define in which ICES statistical rectangle the datapoint is
  ICESrec <- paste(paste(lonmin,lonmax,sep="_"),paste(latmin,latmax,sep="_"), sep="&")
  ## identify positive and negative findings at each ICESrec
  site.occurance <- cbind(ICESrec,sapply(1:length(ICESrec),function(t)
    length(which(ICESrec == ICESrec[t] & my.data$occurrenceStatus == "present"))
  )
  )
  pos.sites <- site.occurance[which(!site.occurance[,2] == "0"),1]
  pos.sites <- unique(pos.sites)
  neg.sites <- site.occurance[which(site.occurance[,2] == "0"),1]
  neg.sites <- unique(neg.sites)
  
  all.site.iters.pos<- list()
  all.site.iters.neg<- list()
  all.site.iters.tot<- list()
  
  set.seed(1)
  #split dataset in five. Repeat the process in five iterations (to check repeatability)
  n.repeat <- 5
  CV.level <- 5
  
  # The number of positive and negative sites (ICES rectangles) are shuffled independently to make sure that there are presences and absences in all datasets
  for(r in 1:n.repeat){
    shuffle.pos <- sample(pos.sites, length(pos.sites), replace = FALSE, prob = NULL)
    shuffle.neg <- sample(neg.sites, length(neg.sites), replace = FALSE, prob = NULL)
    all.site.iters.pos[[r]] <- split(shuffle.pos,  rep(seq(1:CV.level), ceiling(length(pos.sites)/CV.level) )[1:length(pos.sites)]  )
    all.site.iters.neg[[r]] <- split(shuffle.neg,  rep(seq(1:CV.level), ceiling(length(neg.sites)/CV.level) )[1:length(neg.sites)]  )
  }
  all.site.iters.tot <- try(lapply(1:n.repeat,function(r)
    lapply(1:CV.level,function(i)
      union(all.site.iters.pos[[r]][[i]],all.site.iters.neg[[r]][[i]])
    )
  ),silent = TRUE)
  ## get the corresponding occurances
  #length(unlist(all.site.iters.tot[[1]]))
  all.occurance.iters <- try(lapply(1:n.repeat,function(r)
    lapply(1:CV.level,function(i)
      unlist( sapply(1:length(all.site.iters.tot[[r]][[i]]),function(s)
        unlist(my.data$ID[grep(
          all.site.iters.tot[[r]][[i]][s], ICESrec
        )
        ]
        )))
    )
  ),silent = TRUE)
  
  save("all.site.iters.tot",all.site.iters.tot, file= paste(iterations.path,"/site.iters_",species,".rda",sep=""))
  if(length(grep("Error",all.occurance.iters ))==0){
    save("all.occurance.iters", all.occurance.iters, file= paste(iterations.path,"/occurance.iters_",species,".rda",sep=""))
  }
  return(c(length(pos.sites),length(neg.sites)))
}
############################################################
################### calculate ROC curve based on cross validation  ############
##################################################################


calc.ROC <- function(.RF.result,.true.class){
  #..    RF.result[[rep]][["predict.allvar"]]  
  # .RF.result <-  rf.output.cv[[curr.rep]]
   #.true.class <-.true.class
  #.method <-method
  all.allvar.pred <- c()
  
#  if(.method == "allvars"){  
    p.present <- unlist(sapply(1:length(.RF.result),function(i)
      .RF.result[[i]]$probs.sel[,"present"]
    ))
 #   print("all vars")
  #}else{  
#    p.present <- unlist(sapply(1:length(.RF.result),function(i)
  #    .RF.result[[i]]$probs.sel[,"present"]
   # print("sel vars")
 # }
  sorted.values <- .true.class[ as.numeric(names(p.present))]
  #plot(sorted.values+runif(length(sorted.values),-0.2,0.2),p.present,pch=".")
  cutoffs <- seq(1:999)/1000
  sens <- sapply(cutoffs,function(c)
    length(which( p.present[which(sorted.values=="present")]>c ))/
      length(which(sorted.values=="present"))
  )
  spec <-  sapply(cutoffs,function(c)
    length(which( p.present[which(sorted.values=="absent")]<=c ))/
      length(which(sorted.values=="absent"))
  )
  sens <- sens[rev(order(spec))]
  spec <- spec[rev(order(spec))]
  ########### calculate AUC #############
  sens <- c(0,sens,1) #make sure it goes from 0 to 1
  spec <- c(1,spec,0)
  #plot(1-spec,sens)
  one.minus.spec <- 1-spec
  AUC <- 0
  for(i in 2:length(sens)){
    AUC <- AUC + (one.minus.spec[i] - one.minus.spec[i-1])*((sens[i]+sens[i-1])/2)
  }
  # plot((1-spec),sens, pch=".")
  #  lines((1-spec),sens)
  # text(0.5,0.5,AUC)
  return("ROC"=list("sens" = sens,"spec"=spec, "AUC"=AUC))
}

############################################################
################### plot ROC curve   ############
##################################################################

plot.ROC <- function(.ROC.path,.my.species,.mean.ROC,.all.ROC){
  plotname <- paste(ROC.path,"/","plotROC_",.my.species,method,".png",sep="")
  png(plotname)  #to make file
  par(mar = rep(2, 4))

  plot(c(0,0),c(0,-1),xlim=c(0,1),ylim=c(0,1))
##polygon(c((1-mean.spec),1,1),c(mean.sense,1,0),lty=2,lwd=0.2,density=20)
  polygon(c((1-mean.ROC$spec),1,1),c(mean.ROC$sens,1,0),lty=2,lwd=0.2,col=3,border=NA)

  lines(c((1-mean.ROC$spec),1,1,0),c(mean.ROC$sens,1,0,0),lty=2,lwd=3)

  lines(1-mean.ROC$spec,mean.ROC$sens,lty=2,lwd=3)
  #plot(1-spec,sens,xlim=c(0,1),ylim=c(0,1))
  lines(c(0,1),c(0,1))
  # text(0.5,0.95,strsplit(speciesfile,"/")[[1]][2])
  text(0.25,0.9,paste("Mean AUC=",round(mean.ROC$AUC,3)),cex=0.75)
  text(0.5,0.1,"AUC from repetitions",cex=0.75)
  for(r in all.rep){
    lines(1-all.ROC[[r]]$spec,all.ROC[[r]]$sens,lty=2,col=r+1)
    text(0.05+(r/12),0.05,paste(round(all.ROC[[r]]$AUC,3)),cex=0.65,col=1)
  }
  dev.off() 
  return(sapply(1:length(all.ROC), function(r)round(all.ROC[[r]]$AUC,3)))
}
############################################################
################### MCMC.process   ############
##################################################################

MCMC.process = function(.mydata, .iters, .process.plan,.process.ID){
  iter.ID <- .iters[[.process.plan[.process.ID,2]
  ]][[.process.plan[.process.ID,3]
  ]]
 iter.nr <-  unlist(sapply(1:length(iter.ID),function(i)
grep(paste("A",iter.ID[i],"A"),paste("A",.mydata$ID,"A"))
 )) 
 
# length(iter.ID)
# length(iter.nr)
  train <- seq(1:length(.mydata$occurrenceStatus))[-iter.nr ]
  #ncol <- length(names(mydata))# value is the last column
  # temp <- try(randomForest(.mydata[train,-c(ncol)], mydata[train,ncol], prox=TRUE))
  
  # make sure that class is the ast column
  class <- .mydata$occurrenceStatus
  remove.col <- unlist( sapply(c("ID" ,"Lat", "Lon", "occurrenceStatus" ),function(i)
    grep(i,names(.mydata))
  ))
  .mydata$class <- class
 # length(class)
  #unique(class)
  d <- .mydata[train,-remove.col]
  #dim(d)
  #length(train)
  for(i in 1:length(names(d))){
    print(names(d)[i])
    print(unique(d[,i]))
  }
  attribute.names <- names(d)
  names(d) <- c( paste("nr",seq(1,length(d[1,])-1 ),sep="_"),"class")
  d$class <- factor(d$class)
  n.attributes <- length(names(d))-1 
  
  temp <- try(MCMFresult <- mcfs(class~., d, mcfs.projections=600, mcfs.projectionSize=min(2,n.attributes), mcfs.cutoffPermutations=20, mcfs.threadsNumber=8)
  )
  cutoff <- temp$cutoff_value
  RI <- temp$RI
  RI$attribute.name <- attribute.names[ sapply(1:length(RI$attribute),function(i)
    as.numeric( strsplit(RI$attribute[i],"_")[[1]][2]))]
  
  selected.nr <- as.numeric( rownames(temp$RI)[seq(1,temp$cutoff_value)])
 if(length(selected.nr)>0){ 
   selected <- attribute.names[selected.nr]
 }else{
     selected <- NA}
  return(list("cutoff" = cutoff, "RI"=RI,"selected"=selected,"selected.nr"=selected.nr))
}