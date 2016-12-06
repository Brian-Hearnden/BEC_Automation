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
setwd("C:/Users/elilles/Documents/Automated BEC Mapping of Woodland/WoodlandMapping")


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
