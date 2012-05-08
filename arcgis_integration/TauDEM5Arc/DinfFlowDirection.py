# Script Name: DinfFlowDirection
# 
# Created By:  David Tarboton
# Date:        9/23/11

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
ang = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nOutput Dinf Flow Direction File: "+ang)
slp = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nOutput Dinf Slope File: "+slp)

# Construct command 
cmd = 'mpiexec -n ' + inputProc + ' DinfFlowDir -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + ' -slp ' + '"' + slp + '"'
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
arcpy.CalculateStatistics_management(ang)
arcpy.CalculateStatistics_management(slp)
