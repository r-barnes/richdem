# Script Name: DinfDecayingAccumulation
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

decaymultiplier=arcpy.GetParameterAsText(1)
desc = arcpy.Describe(decaymultiplier)
dm=str(desc.catalogPath)
arcpy.AddMessage("\nInput Decay Multiplier Grid: "+dm)

weightgrid=arcpy.GetParameterAsText(2)
if arcpy.Exists(weightgrid):
    desc = arcpy.Describe(weightgrid)
    wg=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Weight Grid: "+wg)

shapefile=arcpy.GetParameterAsText(3)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

edgecontamination=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nEdge Contamination: "+edgecontamination)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
dsca = arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nOutput Decayed Specific Catchment Area Grid: "+dsca)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' DinfDecayAccum -ang ' + '"' + ang + '"' + ' -dsca ' + '"' + dsca + '"' + ' -dm ' + '"' + dm + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if arcpy.Exists(weightgrid):
    cmd = cmd + ' -wg ' + '"' + wg + '"'
if edgecontamination == 'false':
    cmd = cmd + ' -nc '

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
arcpy.CalculateStatistics_management(dsca)
