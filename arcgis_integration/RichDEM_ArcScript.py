# Script Name: Remove Pits
# 
# Created By:  David Tarboton
# Date:        9/21/11

# Import ArcPy site-package and os modules
#
import arcpy
import os
import sys
import time
import string
import subprocess

# Get and describe the first argument
#
inLyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inLyr)
inZfile=str(desc.catalogPath)
arcpy.AddMessage("\nInput Elevation file: "+inZfile)

# Get the Input No. of Processes
#
inputProc=arcpy.GetParameterAsText(1)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Get the output file
#
outFile = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nOutput Pit Removed Elevation file: "+outFile)

# Construct the taudem command line.  Put quotes around file names in case there are spaces
cmd = 'mpiexec -n ' + inputProc + ' pitremove -z ' + '"' + inZfile + '"' + ' -fel ' + '"' + outFile + '"'
arcpy.AddMessage("\nCommand Line: "+cmd)
#os.system(cmd)
#process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
#for line in process.stdout.readlines():
#    arcpy.AddMessage(line)

#  Calculate statistics so that grids display with correct bounds
#arcpy.AddMessage('Executing: Calculate Statistics\n')
#arcpy.CalculateStatistics_management(outFile)
