##### Template R code for mapping of woodland units from known points. 
##### Colin Mahony, UBC Forestry, 778-288-4008, c_mahony@alumni.ubc.ca
##### November 4th, 2016
rm(list=ls()) #clean the workspace so all previous objects are deleted

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkgs = c("scales","MASS", "stats", "rgl", "RColorBrewer", "FNN", "igraph", "raster", "maps"
         , "maptools", "sp", "colorRamps", "rgeos", "rgdal", "foreign", "randomForest")

ipak(pkgs)

## need to create this folder and a "Results" and "InputData" folder in it. 
setwd("C:/Users/elilles/Documents/Automated BEC Mapping of Woodland/WoodlandMapping")

### generate a data frame of analysis variables for the grid points. 
grid.ref <- read.csv("InputData\\ESSFmc_BAFApts4ClimateBC_Normal_1961_1990MSY.csv", strip.white = TRUE, na.strings = c("NA","",-9999) )
nonCNA <- which(is.na(grid.ref[,6]))  # dem cells outside climateNA extent. need this for later

#######WILLS CODE FOR CLIMATE BC VARIABLES TO USE####
##Expects data with PlotNo, BGC, Lat, long, Elev as first 5 columns and ALL Variable output from ClimateWNA
####modify
colnames(X1)[1]=c("BGC")
colnames(X1)[2]=c("_")
records <- nrow(X1)
##X1$PlotNo <- as.integer(seq(from = 1, to = records, by =1))
attr(X1, "row.names") <- (X1$PlotNo)
X2 <- X1 [, c("PlotNo", "BGC", "Latitude", "Longitude", "Elevation")]
X1=X1[,-c(1,3:5)]
# Drop
X1$BGC <- as.factor(X1$BGC)
save(X1,file=paste(fname,".Rda",sep=""))

X1$CMDMax <- X1$CMD07
X1$PPTJune <- X1$PPT06
X1$CMDJune <- X1$CMD06
X1$TMaxJune <- X1$TMax06
X1$DD_18June <- X1$DD_18_06
X1$DD5May <- X1$DD5_05
X1$CMD.grow <- X1$CMD05 + X1$CMD06 +X1$CMD07 +X1$CMD08 +X1$CMD09
X1$PPT.dormant <- X1$PPT_at + X1$PPT_wt
X1$CMD.def <- 500 - (X1$PPT.dormant)
X1$CMD.def [X1$CMD.def < 0] <- 0
X1$CMD.total <- X1$CMD.def + X1$CMD
X1save = X1
#### Choose Biologically Interpretable Variables Only:

TEMP.list=c("Tmax_sp","Tmax_sm","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","DD_0_sp", "DD_0_at",
            "DD5_sp","DD5_sm","DD5_at","DD_18_sp","DD_18_sm","TMaxJune",
            "MAT","MWMT","MCMT","DD5","EMT","EXT", "DD_18June", "DD5May")
PPT.list=c("PPTJune", "PPT_sp", "PPT_sm", "PPT_at", "PPT_wt" ,"MSP", "MAP","PAS", "PPT.dormant")
OTHER.list=c("CMD_sp","CMD_sm","CMDMax","AHM","SHM","NFFD","bFFP","FFP","CMD", "CMD.grow", "CMDJune", "CMD.def", "CMD.total")
ClimateVar=c(TEMP.list,PPT.list,OTHER.list)
List=c("BGC")
X1save = X1

#############Create final data set based on Options selected above
X1$BGC  <- as.factor(X1$BGC)
X1=X1[,names(X1) %in% c(List,ClimateVar)]

save(X1,file=paste(fname,"_DataSet",".Rda",sep=""))

####################
#########END OF WILLS CODE
###########################
#select predictor variables
predictors <- names(grid.ref)[-grep("id|tude|Elev|MAR|RH|18",names(grid.ref))] #remove selected variables from the variable set (could keep elevation in... perhaps it will be a good predictor)
X.grid.ref <- grid.ref[-nonCNA,which(names(grid.ref)%in%predictors)]  #data frame of analysis variables. removes NA values
sum(is.na(X.grid.ref))  #check if there are any NA values. the answer should be "0". 

##log-transform zero-limited variables
zerolim <- grep("MAP|MSP|PAS|DD|CMD|FF|ref",names(X.grid.ref))
for(i in zerolim){X.grid.ref[which(X.grid.ref[,i]==0),i] <- 1}  #set zero values to one, to facilitate log-transformation
X.grid.ref[,zerolim] <- log(X.grid.ref[,zerolim]) #log-transform 

write.csv(X.grid.ref,"InputData\\X.grid.ref.csv", row.names=FALSE)


############
## Analysis
############ 

BGCv10.pts <- droplevels(BGCv10.pts[-nonCNA,]) #remove NA points and unused factor levels
BGC <- BGCv10.pts$MAP_LABEL

## create three classes: subalpine, parkland, and alpine
class <- rep(NA, length(BGC))
class[BGC%in%c("ESSFmcp", "ESSFmkp", "ESSFmvp", "ESSFunp", "ESSFwvp", "MHunp", "MHwhp")] <- "parkland"
class[grep("BAFA|CMA", BGC)] <- "alpine"
class[is.na(class)] <- "subalpine"
class <- factor(class)

# map the three classes
par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
ColScheme <- c("dodgerblue", "yellow", "black")
X <- dem  #uses dem as a template raster
values(X) <- NA
values(X)[land] <- class
plot(X, col=ColScheme, xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
plot(X, col=ColScheme, xaxt="n", yaxt="n", xlim=c(-132, -127.5), ylim=c(56,58), legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 


##################
#### Parkland classification and mapping Trial 1: balanced data
##################

# trial <- "FirstApprox"
# 
# #create a fake set of "known points" by subsampling the grid. in reality, these points will be provided by Will and Erica (and will also require a separte climateNA file)
# training <- sample(1:length(class),10000)
# table(class[training]) 

#classify zone based on plant community

# rf <- randomForest(X.grid.ref[training,], class[training], strata=class[training], sampsize=rep(min(table(class[training])), length(levels(class[training]))))  #train the RF model. the strata and sampsize arguments are used for tree-level downsampling to balance the training set
# rf.pred <- predict(rf, X.grid.ref)  #predict back to the whole grid. 
# ct <- table(group=class,class=rf.pred)
# ClassCorrect <- diag(prop.table(ct, 1))
# AllCorrect <- sum(diag(ct))/sum(ct)
# 
# png(filename=paste("Results\\WoodlandPrediction_",trial,".png",sep=""), type="cairo", units="in", width=12, height=6, pointsize=12, res=400)
# par(mfrow=c(1,2))
# par(mar=c(0,0,1.5,1))
# ColScheme <- c("dodgerblue", "yellow", "black")
# 
# #BGCv10 map
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

#inset zoom map
# xlim=c(-131, -128)
# ylim=c(56.25,57.76)
# rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
# par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
# plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# 
# par(plt = c(0, 1, 0, 0.93), new = F)
# values(X)[land][-nonCNA] <- rf.pred
# plot(X, col=ColScheme,  xaxt="n", yaxt="n", legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 
# box(col="black", lwd=1.5)
# mtext("Random Forest Prediction", 3, adj=0.5, padj=0, cex=1.1, line=0.2)
# legend(extent(X)[1], 58.2, legend=paste(levels(class), ": ", as.integer((1-ClassCorrect)*100), "% error", sep=""), 
#        title=paste("Total error: ", round((1-AllCorrect)*100, 1), "%", sep="") , 
#        cex=1.1, inset=0.01, bty="n")
# text(extent(X)[1], 56.5, paste("change in parkland area: ", round(100*(sum(rf.pred=="parkland")-sum(class=="parkland"))/sum(class=="parkland"),0), "%", sep=""), pos=4)

#inset zoom map
# xlim=c(-131, -128)
# ylim=c(56.25,57.76)
# rect(xlim[1],  ylim[1],  xlim[2],  ylim[2],  col=(alpha("lightgrey", 0.3)))
# par(plt = c(0.01, 0.54, 0.01, 0.5), new = TRUE)
# plot(X, xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, col=ColScheme, legend=FALSE, legend.mar=0, maxpixels=ncell(X)) 

# dev.off()
