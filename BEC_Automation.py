#import libraries
import arcpy
import subprocess
import sys
import os

#tool inputs for R
workspace = GetParameterAsText(0) #workspace type
trainingPoints = GetParameterAsText(1) #Example Arg
outputFC = GetParameterAsText(2) #outputFC name

#tool inputs for python
tolerance = GetParameterAsText(3) #Number
minArea = GetParameterAsText(4) #Number

# get script location
pyScript = sys.argv[0]
toolDir = os.path.dirname(pyScript)
rScriptPath = toolDir + "TBD.r" #Name of script TBD

#Subprocess Args
rCMD = "R --slave --vanilla --args"
args = .join([workspace, trainingPoints, outputFC])

#command to call R
cmd = rCMD + args + "<" + rScriptPath

#Execute command
os.system(cmd)

#project output
sr = arcpy.SpatialReference ("TBD") #Projection TBD
arcpy.management.DefineProjection(outputFC, sr)

params = arcpy.gp.GetParameterInfo()
#renderFile = os.path.join(toolDir, "TBD.shp") #output name
params[2].Symbology = outputFC

#create geodatabase and feature classes
db = arcpy.CreateFileGDB_management(workspace, "BEC_Automation.gdb")
arcpy.FeatureClasstoFeatureClass_conversion (outputFC, db, "orig_BEC")

# tolerance and minimum area must still be determined
arcpy.SimplifyPolygon_cartography(outputFC, db + "\simp_BEC", "BEND_SIMPLIFY", TBD, TBD, "RESOLVE_ERRORS", "NO_KEEP")
