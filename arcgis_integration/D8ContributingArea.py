# Script Name: D8ContributingArea
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

shapefile=arcpy.GetParameterAsText(1)
if arcpy.Exists(shapefile):
    desc = arcpy.Describe(shapefile)
    shfl=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shfl)

weightgrid=arcpy.GetParameterAsText(2)
if arcpy.Exists(weightgrid):
    desc = arcpy.Describe(weightgrid)
    wtgr=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Weight Grid: "+wtgr)

edgecontamination=arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nEdge Contamination: "+edgecontamination)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Output
ad8 = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput D8 Contributing Area Grid: "+ad8)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if arcpy.Exists(weightgrid):
    cmd = cmd + ' -wg ' + '"' + wtgr + '"'
if edgecontamination == 'false':
    cmd = cmd + ' -nc '

##if arcpy.Exists(shapefile) and arcpy.Exists(weightgrid) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -o ' + '"' + shapefile + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -nc '
##elif arcpy.Exists(shapefile) and arcpy.Exists(weightgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -o ' + '"' + shapefile + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -wg ' + '"' + weightgrid + '"'
##elif arcpy.Exists(weightgrid) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -nc '
##elif arcpy.Exists(shapefile) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -o ' + '"' + shapefile + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -nc '
##elif arcpy.Exists(shapefile):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -o ' + '"' + shapefile + '"' + ' -ad8 ' + '"' + ad8 + '"'
##elif arcpy.Exists(weightgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -ad8 ' + '"' + ad8 + '"'
##elif edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"' + ' -nc '
##else:
##    cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"'
    
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
arcpy.CalculateStatistics_management(ad8)



