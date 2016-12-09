rm(list=ls()) #clean the workspace so all previous objects are deleted


#Download projects that are not currently downloaded and load packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkgs = c("scales","MASS", "stats", "rgl", "RColorBrewer", "FNN", "igraph", "raster", "maps"
         , "maptools", "sp", "colorRamps", "rgeos", "rgdal", "foreign", "randomForest")

ipak(pkgs)


setwd("S:\\srm\\smt\\Workarea\\BHearnden\\Ministry_Proj\\Research\\BEC_Automation\\Code_Test")

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


# RUN SAME CODE FOR TRAINING DATASET #

Training <- read.csv("InputData\\Training_Pts4CimateBC_Normal_1961_1990MSY.csv", strip.white = TRUE, na.strings = c("NA","",-9999) )
nonCNA <- which(is.na(Training[,6]))
rownames(Training) <- seq(1, length(Training$id1), 1) 

colnames(Training)[1]=c("BGC")
colnames(Training)[2]=c("id2")
TrainingRecords <- nrow(Training)
Training$id2 <- as.integer(seq(from = 1, to = TrainingRecords, by =1))

Training2 <- Training [, c("BGC", "id2", "Latitude", "Longitude", "Elevation")]
Training=Training[,-c(1,3:5)]

Training$CMDMax <- Training$CMD07
Training$PPTJune <- Training$PPT06
Training$CMDJune <- Training$CMD06
Training$TMaxJune <- Training$TMax06
Training$DD_18June <- Training$DD_18_06
Training$DD5May <- Training$DD5_05
Training$CMD.grow <- Training$CMD05 + Training$CMD06 +Training$CMD07 +Training$CMD08 +Training$CMD09
Training$PPT.dormant <- Training$PPT_at + Training$PPT_wt
Training$CMD.def <- 500 - (Training$PPT.dormant)
Training$CMD.def [Training$CMD.def < 0] <- 0
Training$CMD.total <- Training$CMD.def + Training$CMD
Trainingsave = Training
#### Choose Biologically Interpretable Variables Only:

Training_TEMP.list=c("Tmax_sp","Tmax_sm","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","DD_0_sp", "DD_0_at",
            "DD5_sp","DD5_sm","DD5_at","DD_18_sp","DD_18_sm","TMaxJune",
            "MAT","MWMT","MCMT","DD5","EMT","EXT", "DD_18June", "DD5May")
Training_PPT.list=c("PPTJune", "PPT_sp", "PPT_sm", "PPT_at", "PPT_wt" ,"MSP", "MAP","PAS", "PPT.dormant")
Training_OTHER.list=c("CMD_sp","CMD_sm","CMDMax","AHM","SHM","NFFD","bFFP","FFP","CMD", "CMD.grow", "CMDJune", "CMD.def", "CMD.total")
Training_ClimateVar=c(Training_TEMP.list,Training_PPT.list,Training_OTHER.list)
Training_List=c("BGC")
Trainingsave = Training

#############Create final data set based on Options selected above
Training$BGC  <- as.factor(Training$BGC)
Training=Training[,names(Training) %in% c(List,Training_ClimateVar)]
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
# class <- rep(NA, length(BGC))
# class[BGC%in%c("ESSFmcp", "ESSFmkp", "ESSFmvp", "ESSFunp", "ESSFwvp", "MHunp", "MHwhp")] <- "parkland"
# class[grep("BAFA|CMA", BGC)] <- "alpine"
# class[is.na(class)] <- "subalpine"

class <- ifelse(Training2$BGC == "Woodland", "ESSFmcp", ifelse(Training2$BGC == "ESSFmc", "ESSFmc", "BAFA"))
class <- factor(class)
#trial <- "FirstApprox"

#create a fake set of "known points" by subsampling the grid. in reality, these points will be provided by Will and Erica (and will also require a separte climateNA file)
# training <- read.csv("InputData\\Training_Pts4CimateBC_Normal_1961_1990MSY.csv")
# table(class[training]) 

#classify zone based on plant community
rf <- randomForest(x = Training, y=class, strata=class, sampsize=rep(min(table(class)), length(levels(class))))  #train the RF model. the strata and sampsize arguments are used for tree-level downsampling to balance the Training set
rf.pred <- predict(rf, grid.ref)  #predict back to the whole grid. 
ct <- table(group=class, class=rf.pred)
ClassCorrect <- diag(prop.table(ct, 1))
AllCorrect <- sum(diag(ct))/sum(ct)

