# Script Name: LengthAreaStreamSource
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
plen=str(desc.catalogPath)
arcpy.AddMessage("\nInput Length Grid: "+plen)

inlyr1 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr1)
ad8=str(desc.catalogPath)
arcpy.AddMessage("\nInput Contributing Area Grid: "+ad8)

threshold=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nThreshold(M): "+threshold)

exponent=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nExponent(y): "+exponent)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
ss = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput Stream Source Grid: "+ss)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' LengthArea -plen ' + '"' + plen + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -ss ' + '"' + ss + '"' + ' -par ' + threshold + ' ' + exponent

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
arcpy.CalculateStatistics_management(ss)
