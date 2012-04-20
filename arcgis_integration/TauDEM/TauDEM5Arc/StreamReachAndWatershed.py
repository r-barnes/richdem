# Script Name: StreamReachAndWatershed
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
arcpy.AddMessage("\nInput D8 Drainage Area: "+ad8)

inlyr3 = arcpy.GetParameterAsText(3)
desc = arcpy.Describe(inlyr3)
src=str(desc.catalogPath)
arcpy.AddMessage("\nInput Stream Raster Grid: "+src)

shapefile=arcpy.GetParameterAsText(4)
if arcpy.Exists(shapefile):    
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile as Network Nodes: "+shfl)

delineate=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nDelineate Single Watershed: "+delineate)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
ord = arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nOutput Stream Order Grid: "+ord)

tree = arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nOutput Network Connectivity Tree: "+tree)

coord=arcpy.GetParameterAsText(9)
arcpy.AddMessage("\nOutput Network Coordinates: "+coord)

net=arcpy.GetParameterAsText(10)
arcpy.AddMessage("\nOutput Stream Reach Shapefile: "+net)

w=arcpy.GetParameterAsText(11)
arcpy.AddMessage("\nOutput Watershed Grid: "+w)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' StreamNet -fel ' + '"' + fel + '"' + ' -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -src ' + '"' + src + '"' + ' -ord ' + '"' + ord + '"' + ' -tree ' + '"' + tree + '"' + ' -coord ' + '"' + coord + '"' + ' -net ' + '"' + net + '"' + ' -w ' + '"' + w + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if delineate == 'true':    
    cmd = cmd + ' -sw '

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
arcpy.CalculateStatistics_management(ord)
arcpy.CalculateStatistics_management(net)
arcpy.CalculateStatistics_management(w)
