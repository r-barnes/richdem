# Script Name: D8ExtremeUpslope
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
p=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Flow Direction Grid: "+p)

inlyr2 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr2)
sa=str(desc.catalogPath)
arcpy.AddMessage("\nInput Value Grid: "+sa)

maximumupslope=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nMaximum Upslope: "+maximumupslope)

edgecontamination=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nEdge Contamination: "+edgecontamination)

shapefile=arcpy.GetParameterAsText(4)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
ssa = arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nOutput Extreme Value Grid: "+ssa)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' D8FlowPathExtremeUp -p ' + '"' + p + '"' + ' -sa ' + '"' + sa + '"' + ' -ssa ' + '"' + ssa + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if maximumupslope == 'false':
    cmd = cmd + ' -min '
if edgecontamination == 'false':
    cmd = cmd + ' -nc '

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
arcpy.CalculateStatistics_management(ssa)
