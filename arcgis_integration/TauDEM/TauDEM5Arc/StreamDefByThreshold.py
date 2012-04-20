# Script Name: StreamDefByThreshold
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
ssa=str(desc.catalogPath)
arcpy.AddMessage("\nInput Accumulated Stream Source Grid: "+ssa)

maskgrid = arcpy.GetParameterAsText(1)
if arcpy.Exists(maskgrid):        
    desc = arcpy.Describe(maskgrid)
    mask=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Mask Grid: "+mask)

threshold=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nThreshold: "+threshold)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
src = arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nOutput Stream Raster Grid: "+src)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' Threshold -ssa ' + '"' + ssa + '"' + ' -src ' + '"' + src + '"' + ' -thresh ' + threshold
if arcpy.Exists(maskgrid):
    cmd = cmd + ' -mask ' + mask

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
arcpy.CalculateStatistics_management(src)
