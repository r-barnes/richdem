# Script Name: SlopeAreaCombination
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
slp=str(desc.catalogPath)
arcpy.AddMessage("\nInput Slope Grid: "+slp)

inlyr = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr)
sca=str(desc.catalogPath)
arcpy.AddMessage("\nInput Area Grid: "+sca)

slopeexponent=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nSlope Exponent(m): "+slopeexponent)

areaexponent=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nArea Exponent(n): "+areaexponent)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
sa = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput Slope Area Grid: "+sa)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' SlopeArea -slp ' + '"' + slp + '"' + ' -sca ' + '"' + sca + '"' + ' -sa ' + '"' + sa + '"' + ' -par ' + slopeexponent + ' ' + areaexponent

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
arcpy.CalculateStatistics_management(sa)
