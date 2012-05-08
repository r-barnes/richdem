# Script Name: DinfReverseAccumulation
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
ang=str(desc.catalogPath)
arcpy.AddMessage("\nInput D-Infinity Flow Direction Grid: "+ang)

inlyr1 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr1)
dm=str(desc.catalogPath)
arcpy.AddMessage("\nInput Weight Grid: "+dm)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
racc = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nOutput Reverse Accumulation Grid: "+racc)

dmax = arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nOutput Maximum Downslope Grid: "+dmax)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' DinfRevAccum -ang ' + '"' + ang + '"' + ' -wg ' + '"' + dm + '"' + ' -racc ' + '"' + racc + '"' + ' -dmax ' + '"' + dmax + '"'

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
arcpy.CalculateStatistics_management(racc)
arcpy.CalculateStatistics_management(dmax)
