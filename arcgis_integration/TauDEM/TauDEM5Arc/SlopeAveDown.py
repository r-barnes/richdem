# Script Name: SlopeAveDown
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
fel=str(desc.catalogPath)
arcpy.AddMessage("\nInput Pit Filled Elevation Grid: "+fel)

distance = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nDistance: "+distance)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
slpd = arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nOutput Slope Average Down Grid: "+slpd)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' SlopeAveDown -p ' + '"' + p + '"' + ' -fel ' + '"' + fel + '"' + ' -slpd ' + '"' + slpd + '"' + ' -dn ' + distance

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
arcpy.CalculateStatistics_management(slpd)
