# Script Name: GridNetwork
# 
# Created By:  David Tarboton
# Date:        9/28/11

# Import ArcPy site-package and os modules
import arcpy 
import os
import subprocess

# Inputs
inlyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inlyr)
p=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Flow Direction file: "+p)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(1)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

shapefile=arcpy.GetParameterAsText(2)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

maskgrid=arcpy.GetParameterAsText(3)
if arcpy.Exists(maskgrid):
    desc = arcpy.Describe(maskgrid)
    mkgr=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Mask Grid: "+mkgr)

maskthreshold=arcpy.GetParameterAsText(4)
if maskthreshold:
    arcpy.AddMessage("\nInput Mask Threshold Value: "+maskthreshold)

# Outputs
gord=arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput Strahler Network Order Grid: "+gord)

plen=arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nOutput Longest Upslope Length Grid: "+plen)

tlen = arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nOutput Total Upslope Length Grid: "+tlen)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' GridNet -p ' + '"' + p + '"' + ' -plen ' + '"' + plen + '"' + ' -tlen ' + '"' + tlen + '"' + ' -gord ' + '"' + gord + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if arcpy.Exists(maskgrid):
    cmd = cmd + ' -mask ' + '"' + mkgr + '"' + ' -thresh ' + maskthreshold

##if arcpy.Exists(shapefile) and arcpy.Exists(maskgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' GridNet -p ' + '"' + p + '"' + ' -plen ' + '"' + plen + '"' + ' -tlen ' + '"' + tlen + '"' + ' -gord ' + '"' + gord + '"' + ' -o ' + '"' + shapefile + '"' + ' -mask ' + '"' + maskgrid + '"' + ' -thresh ' + maskthreshold
##elif arcpy.Exists(shapefile):
##    cmd = 'mpiexec -n ' + inputProc + ' GridNet -p ' + '"' + p + '"' + ' -plen ' + '"' + plen + '"' + ' -tlen ' + '"' + tlen + '"' + ' -gord ' + '"' + gord + '"' + ' -o ' + '"' + shapefile + '"'
##elif arcpy.Exists(maskgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' GridNet -p ' + '"' + p + '"' + ' -plen ' + '"' + plen + '"' + ' -tlen ' + '"' + tlen + '"' + ' -gord ' + '"' + gord + '"' + ' -mask ' + '"' + maskgrid + '"' + ' -thresh ' + maskthreshold
##else:
##    cmd = 'mpiexec -n ' + inputProc + ' GridNet -p ' + '"' + p + '"' + ' -plen ' + '"' + plen + '"' + ' -tlen ' + '"' + tlen + '"' + ' -gord ' + '"' + gord + '"'

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
arcpy.CalculateStatistics_management(gord)
arcpy.CalculateStatistics_management(plen)
arcpy.CalculateStatistics_management(tlen)
