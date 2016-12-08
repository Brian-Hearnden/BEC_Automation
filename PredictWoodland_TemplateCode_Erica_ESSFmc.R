##### Template R code for mapping of woodland units from known points. 
##### Colin Mahony, UBC Forestry, 778-288-4008, c_mahony@alumni.ubc.ca
##### November 4th, 2016
rm(list=ls()) #clean the workspace so all previous objects are deleted

library(scales)
library(MASS)   
library(stats)
library(rgl)
library(RColorBrewer)
library(FNN)
library(igraph)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(colorRamps)
library(rgeos)
library(rgdal)
library(foreign)
library(randomForest)

## need to create this folder and a "Results" and "InputData" folder in it. 
setwd("S:\\srm\\smt\\Workarea\\BHearnden\\Ministry_Proj\\Research\\BEC_Automation\\Code_Test")


############
## Data preparation
############ 

## define projections for use in the analysis. 
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.AEA <- CRS ("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs") # standard albers projection for BC gov't data

### BGC v10 linework
ESSFmc_BAFAv10 <- readShapePoly("BEC_SpatialData\\ESSFmc_BAFA_SouthSkeena.shp")
projection(ESSFmc_BAFAv10) <- P4S.AEA # tell R the projection of the shapefile (this is known via the shapefile metadata)
plot(ESSFmc_BAFAv10)
        #Erica added this below to get bec polygons as latlon fo future merge wtih dem pts
        #BGCv10.latlon <- spTransform(BGCv10, P4S.latlon) # reproject to lat-long 


        ### create a study area for the analysis (now the BEC data is only downloaded that contains the study area so this can be skipped)

        ####Skeena region study area
        #regions <- readShapePoly("BEC_SpatialData\\ADM_NR_REG_polygon.shp") # shapefile of BC forest regions
        #projection(regions) <-  P4S.AEA  # tell R the projection of the shapefile (this is known via the shapefile metadata)
        #plot(regions)
        #studyarea <- regions[regions$REGION_NAM=="Skeena Natural Resource Region",]
        #plot(studyarea)
        #studyarea.latlon <- spTransform(studyarea, P4S.latlon) # reproject to lat-long
###ESSF study area
studyarea.latlon<-spTransform(ESSFmc_BAFAv10, P4S.latlon) # reproject to lat-long

### import 30-arcsec digital elevation model (DEM), clip to study area, extract BGCv10 attributes to these points, and create a ClimateNA input file  
dem <- raster("BEC_SpatialData\\namer_dem1.bil")
projection(dem)
dem <- crop(dem, extent(studyarea.latlon)) #crop to study area bounding box (speeds up masking process)
dem <- mask(dem, studyarea.latlon)
land <- which(!is.na(values(dem)))  #use this later for mapping
plot(dem)
pts <- as.data.frame(rasterToPoints(dem)) #convert raster to points
coordinates(pts) <- pts[,1:2]   #promote data frame to spatial points object
projection(pts) <- P4S.latlon  #redefine as latlon

pts2 <- spTransform(pts, P4S.AEA) #project to albers
ESSFmc_BAFAv10.pts <- over(pts2,ESSFmc_BAFAv10)  ##extract the BGCv10 attributes to the points (this takes a while, like 30ish minutes)
 write.csv(ESSFmc_BAFAv10.pts,"InputData\\BGCv10.pts.csv", row.names=FALSE) #write to csv so that you don't have to do the overlay more than once. 

ESSFmc_BAFAv10.pts <- read.csv("InputData\\BGCv10.pts.csv") #if this file has not already been created, then make it with the two rows above. 
pts.data <- as.data.frame(pts2)
pts.data1 <- data.frame(id1=ESSFmc_BAFAv10.pts$MAP_LABEL, id2=ESSFmc_BAFAv10.pts$PHASE, lat=pts.data$y, lon=pts.data$x, el=pts.data$namer_dem1) ## create the climateNA input file
write.csv(pts.data1,"InputData\\ESSFmc_BAFApts4ClimateBC.csv", row.names=FALSE)

#####NOTE: on file version- remove intermediate file exports (write.csv)

### NON-CODED PROCESS: generate climate data for the pts.csv file using climate NA mult method (annual variables, 1971-2000 normals). 

### generate a data frame of analysis variables for the grid points. 
grid.ref <- read.csv("InputData\\ESSFmc_BAFApts4ClimateBC_Normal_1961_1990MSY.csv", strip.white = TRUE, na.strings = c("NA","",-9999) )
nonCNA <- which(is.na(grid.ref[,6]))
rownames(grid.ref) <- seq(1, length(grid.ref$id1), 1) 

# dem cells outside climateNA extent. need this for later

#######WILLS CODE FOR CLIMATE BC VARIABLES TO USE####
##Expects data with UnqiueID, BGC, Lat, long, Elev as first 5 columns and ALL Variable output from ClimateWNA
####modify
colnames(grid.ref)[1]=c("BGC")
colnames(grid.ref)[2]=c("id2")
records <- nrow(grid.ref)
grid.ref$id2 <- as.integer(seq(from = 1, to = records, by =1))
#attr(grid.ref, "row.names") <- (grid.ref$UnqiueID)
#rownames(grid.ref) <- grid.ref$UniqueID

grid.ref2 <- grid.ref [, c("BGC", "id2", "Latitude", "Longitude", "Elevation")]
grid.ref=grid.ref[,-c(1,3:5)]
# Drop
#grid.ref$BGC <- as.factor(grid.ref$BGC)
#save(grid.ref,file=paste(fname,".Rda",sep=""))

grid.ref$CMDMax <- grid.ref$CMD07
grid.ref$PPTJune <- grid.ref$PPT06
grid.ref$CMDJune <- grid.ref$CMD06
grid.ref$TMaxJune <- grid.ref$TMax06
grid.ref$DD_18June <- grid.ref$DD_18_06
grid.ref$DD5May <- grid.ref$DD5_05
grid.ref$CMD.grow <- grid.ref$CMD05 + grid.ref$CMD06 +grid.ref$CMD07 +grid.ref$CMD08 +grid.ref$CMD09
grid.ref$PPT.dormant <- grid.ref$PPT_at + grid.ref$PPT_wt
grid.ref$CMD.def <- 500 - (grid.ref$PPT.dormant)
grid.ref$CMD.def [grid.ref$CMD.def < 0] <- 0
grid.ref$CMD.total <- grid.ref$CMD.def + grid.ref$CMD
grid.refsave = grid.ref
#### Choose Biologically Interpretable Variables Only:

TEMP.list=c("Tmax_sp","Tmax_sm","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","DD_0_sp", "DD_0_at",
            "DD5_sp","DD5_sm","DD5_at","DD_18_sp","DD_18_sm","TMaxJune",
            "MAT","MWMT","MCMT","DD5","EMT","EXT", "DD_18June", "DD5May")
PPT.list=c("PPTJune", "PPT_sp", "PPT_sm", "PPT_at", "PPT_wt" ,"MSP", "MAP","PAS", "PPT.dormant")
OTHER.list=c("CMD_sp","CMD_sm","CMDMax","AHM","SHM","NFFD","bFFP","FFP","CMD", "CMD.grow", "CMDJune", "CMD.def", "CMD.total")
ClimateVar=c(TEMP.list,PPT.list,OTHER.list)
List=c("BGC")
grid.refsave = grid.ref

#############Create final data set based on Options selected above
grid.ref$BGC  <- as.factor(grid.ref$BGC)
grid.ref=grid.ref[,names(grid.ref) %in% c(List,ClimateVar)]

save(grid.ref,file=paste(fname,"_DataSet",".Rda",sep=""))

####################
#########END OF WILLS CODE
###########################
#select predictor variables
# predictors <- names(grid.ref)[-grep("id|tude|Elev|MAR|RH|18",names(grid.ref))] #remove selected variables from the variable set (could keep elevation in... perhaps it will be a good predictor)
# X.grid.ref <- grid.ref[-nonCNA,which(names(grid.ref)%in%predictors)]  #data frame of analysis variables. removes NA values
# sum(is.na(X.grid.ref))  #check if there are any NA values. the answer should be "0". 

# ##log-transform zero-limited variables
# zerolim <- grep("MAP|MSP|PAS|DD|CMD|FF|ref",names(X.grid.ref))
# for(i in zerolim){X.grid.ref[which(X.grid.ref[,i]==0),i] <- 1}  #set zero values to one, to facilitate log-transformation
# X.grid.ref[,zerolim] <- log(X.grid.ref[,zerolim]) #log-transform 
# 
# write.csv(X.grid.ref,"InputData\\X.grid.ref.csv", row.names=FALSE)


############
## Analysis
############ 

# BGCv10.pts <- droplevels(BGCv10.pts[-nonCNA,]) #remove NA points and unused factor levels
BGC <- grid.ref2$BGC
# 
# ## create three classes: subalpine, parkland, and alpine
class <- rep(NA, length(BGC))
class[BGC%in%c("ESSFmcp", "ESSFmkp", "ESSFmvp", "ESSFunp", "ESSFwvp", "MHunp", "MHwhp")] <- "parkland"
class[grep("BAFA|CMA", BGC)] <- "alpine"
class[is.na(class)] <- "subalpine"
class <- factor(class)
# 
# # map the three classes
# par(mfrow=c(1,1))
# par(mar=c(0,0,0,0))
# ColScheme <- c("dodgerblue", "yellow", "black")
# X <- dem  #uses dem as a template raster
# values(X) <- NA
# values(X)[land] <- class
# plot(X, col=ColScheme, xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# plot(X, col=ColScheme, xaxt="n", yaxt="n", xlim=c(-132, -127.5), ylim=c(56,58), legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# 

##################
#### Parkland classification and mapping Trial 1: balanced data
##################

trial <- "FirstApprox"

#create a fake set of "known points" by subsampling the grid. in reality, these points will be provided by Will and Erica (and will also require a separte climateNA file)
training <- sample(1:length(class),10000)
table(class[training]) 

#classify zone based on plant community
rf <- randomForest(grid.ref[training,], class[training], strata=class[training], sampsize=rep(min(table(class[training])), length(levels(class[training]))))  #train the RF model. the strata and sampsize arguments are used for tree-level downsampling to balance the training set
rf.pred <- predict(rf, grid.ref)  #predict back to the whole grid. 
ct <- table(group=class,class=rf.pred)
ClassCorrect <- diag(prop.table(ct, 1))
AllCorrect <- sum(diag(ct))/sum(ct)

# png(filename=paste("Results\\WoodlandPrediction_",trial,".png",sep=""), type="cairo", units="in", width=12, height=6, pointsize=12, res=400)
# par(mfrow=c(1,2))
# par(mar=c(0,0,1.5,1))
# ColScheme <- c("dodgerblue", "yellow", "black")

#BGCv10 map
# par(plt = c(0, 1, 0, 0.93), new = F)
# X <- dem
# values(X) <- NA
# values(X)[land][-nonCNA] <- class
# plot(X, col=ColScheme, xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# box(col="black", lwd=1.5)
# mtext("BGC mapping v10", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
# legend(extent(X)[1]-0.25, 58, legend=paste(levels(class), ": n=", table(class[training]), " (", round(100*table(class[training])/table(class),0), "%)", sep=""), 
#        title=paste("Training sample: n=", sum(table(class[training])), sep="") , 
#        fill=ColScheme, bg="white", col="lightgrey", box.col="white", cex=1.1, inset=0.01)
# 
# #inset zoom map
# xlim=c(-131, -128)
# ylim=c(56.25,57.76)
# rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
# par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
# plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

par(plt = c(0, 1, 0, 0.93), new = F)
values(X)[land][-nonCNA] <- rf.pred
plot(X, col=ColScheme,  xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
box(col="black", lwd=1.5)
mtext("Random Forest Prediction", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
legend(extent(X)[1], 58.2, legend=paste(levels(class), ": ", as.integer((1-ClassCorrect)*100), "% error", sep=""), 
       title=paste("Total error: ", round((1-AllCorrect)*100, 1), "%", sep="") , 
       cex=1.1, inset=0.01, bty="n")
text(extent(X)[1], 56.5, paste("change in parkland area: ", round(100*(sum(rf.pred=="parkland")-sum(class=="parkland"))/sum(class=="parkland"),0), "%", sep=""), pos=4)

#inset zoom map
xlim=c(-131, -128)
ylim=c(56.25,57.76)
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

dev.off()


##################
#### Parkland classification and mapping Trial 2: unbalanced data
##################

trial <- "Unbalanced"

#create a fake set of "known points" by subsampling the grid. in reality, these points will be provided by Will and Erica (and will also require a separte climateNA file)
training <- sample(1:length(class),10000)
table(class[training]) 

#classify zone based on plant community
rf <- randomForest(X.grid.ref[training,], class[training])  #train the RF model without balancing the training set
rf.pred <- predict(rf, X.grid.ref)  #predict back to the whole grid. 
ct <- table(group=class,class=rf.pred)
ClassCorrect <- diag(prop.table(ct, 1))
AllCorrect <- sum(diag(ct))/sum(ct)

png(filename=paste("Results\\WoodlandPrediction_",trial,".png",sep=""), type="cairo", units="in", width=12, height=6, pointsize=12, res=400)
par(mfrow=c(1,2))
par(mar=c(0,0,1.5,1))
ColScheme <- c("dodgerblue", "yellow", "black")

#BGCv10 map
par(plt = c(0, 1, 0, 0.93), new = F)
X <- dem
values(X) <- NA
values(X)[land][-nonCNA] <- class
plot(X, col=ColScheme, xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
box(col="black", lwd=1.5)
mtext("BGC mapping v10", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
legend(extent(X)[1]-0.25, 58, legend=paste(levels(class), ": n=", table(class[training]), " (", round(100*table(class[training])/table(class),0), "%)", sep=""), 
       title=paste("Training sample: n=", sum(table(class[training])), sep="") , 
       fill=ColScheme, bg="white", col="lightgrey", box.col="white", cex=1.1, inset=0.01)

#inset zoom map
xlim=c(-131, -128)
ylim=c(56.25,57.76)
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

par(plt = c(0, 1, 0, 0.93), new = F)
values(X)[land][-nonCNA] <- rf.pred
plot(X, col=ColScheme,  xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
box(col="black", lwd=1.5)
mtext("Random Forest Prediction", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
legend(extent(X)[1], 58.2, legend=paste(levels(class), ": ", as.integer((1-ClassCorrect)*100), "% error", sep=""), 
       title=paste("Total error: ", round((1-AllCorrect)*100, 1), "%", sep="") , 
       cex=1.1, inset=0.01, bty="n")
text(extent(X)[1], 56.5, paste("change in parkland area: ", round(100*(sum(rf.pred=="parkland")-sum(class=="parkland"))/sum(class=="parkland"),0), "%", sep=""), pos=4)

#inset zoom map
xlim=c(-131, -128)
ylim=c(56.25,57.76)
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

dev.off()




##################
#### Parkland classification and mapping Trial 3: smaller sample size
##################

trial <- "Unbalanced1000"

#create a fake set of "known points" by subsampling the grid. in reality, these points will be provided by Will and Erica (and will also require a separte climateNA file)
training <- sample(1:length(class),1000)
table(class[training]) 

#classify zone based on plant community
rf <- randomForest(X.grid.ref[training,], class[training])  #train the RF model without balancing the training set
rf.pred <- predict(rf, X.grid.ref)  #predict back to the whole grid. 
ct <- table(group=class,class=rf.pred)
ClassCorrect <- diag(prop.table(ct, 1))
AllCorrect <- sum(diag(ct))/sum(ct)

png(filename=paste("Results\\WoodlandPrediction_",trial,".png",sep=""), type="cairo", units="in", width=12, height=6, pointsize=12, res=400)
par(mfrow=c(1,2))
par(mar=c(0,0,1.5,1))
ColScheme <- c("dodgerblue", "yellow", "black")


#BGCv10 map
par(plt = c(0, 1, 0, 0.93), new = F)
X <- dem
values(X) <- NA
values(X)[land][-nonCNA] <- class
plot(X, col=ColScheme, xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
box(col="black", lwd=1.5)
mtext("BGC mapping v10", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
legend(extent(X)[1]-0.25, 58, legend=paste(levels(class), ": n=", table(class[training]), " (", round(100*table(class[training])/table(class),0), "%)", sep=""), 
       title=paste("Training sample: n=", sum(table(class[training])), sep="") , 
       fill=ColScheme, bg="white", col="lightgrey", box.col="white", cex=1.1, inset=0.01)

#inset zoom map
xlim=c(-131, -128)
ylim=c(56.25,57.76)
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

par(plt = c(0, 1, 0, 0.93), new = F)
values(X)[land][-nonCNA] <- rf.pred
plot(X, col=ColScheme,  xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
box(col="black", lwd=1.5)
mtext("Random Forest Prediction", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
legend(extent(X)[1], 58.2, legend=paste(levels(class), ": ", as.integer((1-ClassCorrect)*100), "% error", sep=""), 
       title=paste("Total error: ", round((1-AllCorrect)*100, 1), "%", sep="") , 
       cex=1.1, inset=0.01, bty="n")
text(extent(X)[1], 56.5, paste("change in parkland area: ", round(100*(sum(rf.pred=="parkland")-sum(class=="parkland"))/sum(class=="parkland"),0), "%", sep=""), pos=4)

#inset zoom map
xlim=c(-131, -128)
ylim=c(56.25,57.76)
rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

dev.off()

