# Script Name: D8FlowDirection
# 
# Created By:  David Tarboton
# Date:        9/22/11

# Import ArcPy site-package and os modules
import arcpy 
import os
import subprocess

# Input
inlyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inlyr)
fel=str(desc.catalogPath)
arcpy.AddMessage("\nInput Pit Filled Elevation file: "+fel)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(1)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
p = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nOutput D8 Flow Direction File: "+p)
sd8 = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nOutput D8 Slope File: "+sd8)

# Construct command 
cmd = 'mpiexec -n ' + inputProc + ' D8FlowDir -fel ' + '"' + fel + '"' + ' -p ' + '"' + p + '"' + ' -sd8 ' + '"' + sd8 + '"'
arcpy.AddMessage("\nCommand Line: "+cmd)

# Submit command to operating system
os.system(cmd)

# Capture the contents of shell command and print it to the arcgis dialog box
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)

# Calculate statistics on the output so that it displays properly
arcpy.AddMessage('Executing: Calculate Statistics\n')
arcpy.CalculateStatistics_management(p)
arcpy.CalculateStatistics_management(sd8)
