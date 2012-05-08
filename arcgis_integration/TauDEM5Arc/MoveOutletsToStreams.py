# Script Name: MoveOuletsToStreams
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

inlyr1 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr1)
src=str(desc.catalogPath)
arcpy.AddMessage("\nInput Stream Raster Grid: "+src)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
shfl=str(desc.catalogPath)
arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

maxdistance=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nMinimum Threshold Value: "+maxdistance)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
om = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput Outlet Shapefile: "+om)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' MoveOutletsToStreams -p ' + '"' + p + '"' + ' -src ' + '"' + src + '"' + ' -o ' + '"' + shfl + '"' + ' -om ' + '"' + om + '"' + ' -md ' + maxdistance

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
arcpy.CalculateStatistics_management(om)
