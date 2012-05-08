# Script Name: DinfContributingArea
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
ang=str(desc.catalogPath)
arcpy.AddMessage("\nInput Dinf Flow Direction file: "+ang)

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
sca = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nOutput Dinf Specific Catchment Area Grid: "+sca)

# Construct command
cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -sca ' + '"' + sca + '"'
if arcpy.Exists(shapefile):
    cmd = cmd + ' -o ' + '"' + shfl + '"'
if arcpy.Exists(weightgrid):
    cmd = cmd + ' -wg ' + '"' + wtgr + '"'
if edgecontamination == 'false':
    cmd = cmd + ' -nc '

##if arcpy.Exists(shapefile) and arcpy.Exists(weightgrid) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -o ' + '"' + shapefile + '"' + ' -sca ' + '"' + sca + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -nc '
##elif arcpy.Exists(shapefile) and arcpy.Exists(weightgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -o ' + '"' + shapefile + '"' + ' -sca ' + '"' + sca + '"' + ' -wg ' + '"' + weightgrid + '"'
##elif arcpy.Exists(weightgrid) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -sca ' + '"' + sca + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -nc '
##elif arcpy.Exists(shapefile) and edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -o ' + '"' + shapefile + '"' + ' -sca ' + '"' + sca + '"' + ' -nc '
##elif arcpy.Exists(shapefile):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -o ' + '"' + shapefile + '"' + ' -sca ' + '"' + sca + '"'
##elif arcpy.Exists(weightgrid):
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -wg ' + '"' + weightgrid + '"' + ' -sca ' + '"' + sca + '"'
##elif edgecontamination == 'false':
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -sca ' + '"' + sca + '"' + ' -nc '
##else:
##    cmd = 'mpiexec -n ' + inputProc + ' AreaDinf -ang ' + '"' + ang + '"' + ' -sca ' + '"' + sca + '"'
    
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
arcpy.CalculateStatistics_management(sca)
