# Script Name: DinfTransLimitAccum
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
tsup=str(desc.catalogPath)
arcpy.AddMessage("\nInput Supply Grid: "+tsup)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
tc=str(desc.catalogPath)
arcpy.AddMessage("\nInput Transport Capacity Grid: "+tc)

inlyr3 = arcpy.GetParameterAsText(3)
if arcpy.Exists(inlyr3):
    desc = arcpy.Describe(inlyr3)
    cs=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Concentration Grid: "+cs)

shapefile=arcpy.GetParameterAsText(4)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

edgecontamination=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nEdge Contamination: "+edgecontamination)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
tla = arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nOutput Transport Limited Accumulation Grid: "+tla)

tdep = arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nOutput Deposition Grid: "+tdep)

ctpt = arcpy.GetParameterAsText(9)
if arcpy.Exists(inlyr3):
    arcpy.AddMessage("\nOutput Concentration Grid: "+ctpt)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' DinfTransLimAccum -ang ' + '"' + ang + '"' + ' -tsup ' + '"' + tsup + '"' + ' -tc ' + '"' + tc + '"' + ' -tla ' + '"' + tla + '"' + ' -tdep ' + '"' + tdep + '"'
if arcpy.Exists(inlyr3):
    cmd = cmd + ' -cs ' + '"' + cs + '"' + ' -ctpt ' + '"' + ctpt + '"'
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
arcpy.CalculateStatistics_management(tla)
arcpy.CalculateStatistics_management(tdep)
if arcpy.Exists(inlyr3):
    arcpy.CalculateStatistics_management(ctpt)
