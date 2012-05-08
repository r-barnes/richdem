# Script Name: DinfConcenLimitAccum
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
wg=str(desc.catalogPath)
arcpy.AddMessage("\nInput Effective Runoff Weight Grid: "+wg)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
dg=str(desc.catalogPath)
arcpy.AddMessage("\nInput Disturbance Indicator Grid: "+dg)

inlyr3 = arcpy.GetParameterAsText(3)
desc = arcpy.Describe(inlyr3)
dm=str(desc.catalogPath)
arcpy.AddMessage("\nInput Decay Multiplier Grid: "+dm)

shapefile=arcpy.GetParameterAsText(4)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

concthresh=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nConcentration Threshold: "+concthresh)

edgecontamination=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nEdge Contamination: "+edgecontamination)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
q = arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nOutput Overland Flow Specific Discharge Grid: "+q)

ctpt = arcpy.GetParameterAsText(9)
arcpy.AddMessage("\nOutput Concentration Grid: "+ctpt)

# Construct command 1
cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -sca ' + '"' + q + '"' + ' -wg ' + '"' + wg + '"'
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
arcpy.CalculateStatistics_management(q)

# Construct command 2
cmd = 'mpiexec -n ' + inputProc + ' DinfConcLimAccum -ang ' + '"' + ang + '"' + ' -dg ' + '"' + dg + '"' + ' -dm ' + '"' + dm + '"' + ' -ctpt ' + '"' + ctpt + '"' + ' -q ' + '"' + q + '"' + ' -csol ' + concthresh
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
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
arcpy.CalculateStatistics_management(ctpt)
