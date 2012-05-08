# Script Name: SlopeOverAreaRatio
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
slp=str(desc.catalogPath)
arcpy.AddMessage("\nInput Slope Grid: "+slp)

inlyr2 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr2)
sca=str(desc.catalogPath)
arcpy.AddMessage("\nInput Secific Catchment Area Grid: "+sca)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
sar = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nOutput Slope Divided By Area Ratio Grid: "+sar)

# Construct command 
cmd = 'mpiexec -n ' + inputProc + ' SlopeAreaRatio -slp ' + '"' + slp + '"' + ' -sca ' + '"' + sca + '"' + ' -sar ' + '"' + sar + '"'
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
arcpy.CalculateStatistics_management(sar)
