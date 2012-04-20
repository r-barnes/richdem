# Script Name: PeukarDouglasStreamDef
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
p=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Flow Direction Grid: "+p)

weightcenter = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nWeight Center: "+weightcenter)

weightside = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nWeight Side: "+weightside)

weightdiag = arcpy.GetParameterAsText(4)
arcpy.AddMessage("\nWeight Diagonal: "+weightdiag)

accthresh = arcpy.GetParameterAsText(5)
arcpy.AddMessage("\nAccumulation Threshold: "+accthresh)

contcheck = arcpy.GetParameterAsText(6)
arcpy.AddMessage("\nEdge contamination checking: "+contcheck)

shapelyr=arcpy.GetParameterAsText(7)
if arcpy.Exists(shapelyr):
    desc = arcpy.Describe(shapelyr)
    shapefile=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shapefile)

masklyr=arcpy.GetParameterAsText(8)
if arcpy.Exists(masklyr):
    desc = arcpy.Describe(masklyr)
    mask=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Mask Grid: "+mask)

ad8lyr=arcpy.GetParameterAsText(9)
if arcpy.Exists(ad8lyr):
    desc = arcpy.Describe(ad8lyr)
    ad8=str(desc.catalogPath)
    arcpy.AddMessage("\nInput D8 Contributing Area for Drop Analysis: "+ad8)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(10)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
ss = arcpy.GetParameterAsText(11)
arcpy.AddMessage("\nOutput Stream Source Grid: "+ss)

ssa = arcpy.GetParameterAsText(12)
arcpy.AddMessage("\nOutput Accumulated Stream Source Grid: "+ssa)

src = arcpy.GetParameterAsText(13)
arcpy.AddMessage("\nOutput Stream Raster Grid: "+src)

drp=arcpy.GetParameterAsText(14)
if arcpy.Exists(drp):    
    arcpy.AddMessage("\nOutput Drop Analysis Table: "+drp)

usedroprange = arcpy.GetParameterAsText(15)
arcpy.AddMessage("\nSelect Threshold by Drop Analysis: "+usedroprange)

minthresh=arcpy.GetParameterAsText(16)
arcpy.AddMessage("\nMinimum Threshold Value: "+minthresh)

maxthresh=arcpy.GetParameterAsText(17)
arcpy.AddMessage("\nMaximum Threshold Value: "+maxthresh)

numthresh=arcpy.GetParameterAsText(18)
arcpy.AddMessage("\nNumber of Threshold Values: "+numthresh)

logspace=arcpy.GetParameterAsText(19)
arcpy.AddMessage("\nLogarithmic Spacing: "+logspace)

# Construct first command
cmd = 'mpiexec -n ' + inputProc + ' PeukerDouglas -fel ' + '"' + fel + '"' + ' -ss ' + '"' + ss + '"' + ' -par ' + weightcenter + ' ' + weightside + ' ' + weightdiag 
arcpy.AddMessage("\nCommand Line: "+cmd)
# Submit command to operating system
os.system(cmd)
# Capture the contents of shell command and print it to the arcgis dialog box
process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)
arcpy.CalculateStatistics_management(ss)

# Construct second command
cmd = 'mpiexec -n ' + inputProc + ' AreaD8 -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ssa + '"'
cmd = cmd + ' -wg ' + '"' + ss + '"'
if arcpy.Exists(shapelyr):
    cmd = cmd + ' -o ' + '"' + shapefile + '"'
if contcheck == 'false':
    cmd = cmd + ' -nc'
arcpy.AddMessage("\nCommand Line: "+cmd)
# Submit command to operating system
os.system(cmd)
# Capture the contents of shell command and print it to the arcgis dialog box
process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)
arcpy.CalculateStatistics_management(ssa)

if (usedroprange == 'true') and arcpy.Exists(shapelyr):
    # Construct third command
    cmd = 'mpiexec -n ' + inputProc + ' DropAnalysis -fel ' + '"' + fel + '"' + ' -p ' + '"' + p + '"' + ' -ad8 ' + '"' + ad8 + '"'
    cmd = cmd + ' -ssa ' + '"' + ssa + '"' + ' -o ' + '"' + shapefile + '"' + ' -drp ' + '"' + drp + '"' + ' -par ' + minthresh + ' ' + maxthresh + ' ' + numthresh + ' '
    if logspace == 'false':    
        cmd = cmd + '1'
    else:
        cmd = cmd + '0'
        
    arcpy.AddMessage("\nCommand Line: "+cmd)
    # Submit command to operating system
    os.system(cmd)
    # Capture the contents of shell command and print it to the arcgis dialog box
    process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    arcpy.AddMessage('\nProcess started:\n')
    for line in process.stdout.readlines():
        arcpy.AddMessage(line)
        #(threshold,rest)=line.split(' ',1)

    drpfile = open(drp,"r")
    theContents=drpfile.read()
    (beg,threshold)=theContents.rsplit(' ',1)
    drpfile.close()

# Construct fourth command
cmd = 'mpiexec -n ' + inputProc + ' Threshold -ssa ' + '"' + ssa + '"' + ' -src ' + '"' + src + '"'
if (usedroprange == 'true') and arcpy.Exists(shapelyr):
    cmd = cmd + ' -thresh ' + threshold
else:
    cmd = cmd + ' -thresh ' + accthresh    
if arcpy.Exists(masklyr):
    cmd = cmd + ' -mask ' + '"' + mask + '"'
arcpy.AddMessage("\nCommand Line: "+cmd)
# Submit command to operating system
os.system(cmd)
# Capture the contents of shell command and print it to the arcgis dialog box
process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)
arcpy.CalculateStatistics_management(src)
