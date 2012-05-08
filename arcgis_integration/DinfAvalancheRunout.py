# Script Name: DinfAvalancheRunout
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
ang=str(desc.catalogPath)
arcpy.AddMessage("\nInput D-Infinity Flow Direction Grid: "+ang)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
ass=str(desc.catalogPath)
arcpy.AddMessage("\nInput Avalanche Source Site Grid: "+ass)

propthresh=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nInput Proportion Threshold: "+propthresh)

alphthresh=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nInput Alpha Angle Threshold: "+alphthresh)

pathdistance=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nPath Distance Method: "+pathdistance)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
rz = arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nOutput Runout Zone Grid: "+rz)

dfs = arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nOutput Path Distance Grid: "+dfs)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' DinfAvalanche -fel ' + '"' + fel + '"' + ' -ang ' + '"' + ang + '"' + ' -ass ' + '"' + ass + '"' + ' -rz ' + '"' + rz + '"' + ' -dfs ' + '"' + dfs + '"' + ' -thresh ' + propthresh + ' -alpha ' + alphthresh
if pathdistance == 'Straight Line':
    cmd = cmd + ' -direct '

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
arcpy.CalculateStatistics_management(rz)
arcpy.CalculateStatistics_management(dfs)
