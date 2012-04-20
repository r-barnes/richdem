# Script Name: StreamDropAnalysis
# 
# Created By:  David Tarboton
# Date:        9/29/11

# Import ArcPy site-package and os modules
import arcpy 
import os
import subprocess

# Inputs
inlyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inlyr)
fel=str(desc.catalogPath)
arcpy.AddMessage("\nInput Pit Filled Elevation Grid: "+fel)

inlyr1 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr1)
p=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Flow Direction Grid: "+p)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
ad8=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Contributing Area Grid: "+ad8)

inlyr3 = arcpy.GetParameterAsText(3)
desc = arcpy.Describe(inlyr3)
ssa=str(desc.catalogPath)
arcpy.AddMessage("\nInput Accumulated Stream Source Grid: "+ssa)

shapefile=arcpy.GetParameterAsText(4)
desc = arcpy.Describe(shapefile)
shfl=str(desc.catalogPath)
arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

minthresh=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nMinimum Threshold Value: "+minthresh)

maxthresh=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nMaximum Threshold Value: "+maxthresh)

numthresh=arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nNumber of Threshold Values: "+numthresh)

logspace=arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nLogarithmic Spacing: "+logspace)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(9)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
drp = arcpy.GetParameterAsText(10)
arcpy.AddMessage("\nOutput Drop Analysis Text File: "+drp)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' DropAnalysis -fel ' + '"' + fel + '"' + ' -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -ssa ' + '"' + ssa + '"' + ' -o ' + '"' + shfl + '"' + ' -drp ' + '"' + drp + '"' + ' -par ' + minthresh + ' ' + maxthresh + ' ' + numthresh + ' '
if logspace == 'false':    
    cmd = cmd + '1'
else:
    cmd = cmd + '0'

arcpy.AddMessage("\nCommand Line: "+cmd)

# Submit command to operating system
os.system(cmd)

# Capture the contents of shell command and print it to the arcgis dialog box
process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)

# Calculate statistics on the output so that it displays properly
arcpy.AddMessage('Executing: Calculate Statistics\n')
arcpy.CalculateStatistics_management(drp)
