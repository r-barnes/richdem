# Script Name: SlopeAreaStreamDefinition
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
p=str(desc.catalogPath)
arcpy.AddMessage("\nInput D8 Flow Direction Grid: "+p)

inlyr1 = arcpy.GetParameterAsText(1)
desc = arcpy.Describe(inlyr1)
sca=str(desc.catalogPath)
arcpy.AddMessage("\nInput D-Infinity Contributing Area Grid: "+sca)

inlyr2 = arcpy.GetParameterAsText(2)
desc = arcpy.Describe(inlyr2)
slp=str(desc.catalogPath)
arcpy.AddMessage("\nInput Slope Grid: "+slp)

masklyr=arcpy.GetParameterAsText(3)
if arcpy.Exists(masklyr):
    desc = arcpy.Describe(masklyr)
    mask=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Mask Grid: "+mask)

shapelyr=arcpy.GetParameterAsText(4)
if arcpy.Exists(shapelyr):
    desc = arcpy.Describe(shapelyr)
    shapefile=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Outlets Shapefile: "+shapefile)

fellyr=arcpy.GetParameterAsText(5)
if arcpy.Exists(fellyr):
    desc = arcpy.Describe(fellyr)
    fel=str(desc.catalogPath)
    arcpy.AddMessage("\nInput Pit Filled Elevation Grid for Drop Analysis: "+fel)

ad8lyr=arcpy.GetParameterAsText(6)
if arcpy.Exists(ad8lyr):
    desc = arcpy.Describe(ad8lyr)
    ad8=str(desc.catalogPath)
    arcpy.AddMessage("\nInput D8 Contributing Area Grid for Drop Analysis: "+ad8)

contcheck = arcpy.GetParameterAsText(7)
arcpy.AddMessage("\nEdge contamination checking: "+contcheck)

slpexp = arcpy.GetParameterAsText(8)
arcpy.AddMessage("\nSlope Exponent(m): "+slpexp)

areaexp = arcpy.GetParameterAsText(9)
arcpy.AddMessage("\nArea Exponent(n): "+areaexp)

thresh = arcpy.GetParameterAsText(10)
arcpy.AddMessage("\nThreshold(T): "+thresh)

# Input Number of Processes
inputProc=arcpy.GetParameterAsText(11)
arcpy.AddMessage("\nInput Number of Processes: "+inputProc)

# Outputs
src = arcpy.GetParameterAsText(12)
arcpy.AddMessage("\nOutput Stream Raster Grid: "+src)

sa = arcpy.GetParameterAsText(13)
arcpy.AddMessage("\nOutput Source Area Grid: "+sa)

ssa = arcpy.GetParameterAsText(14)
arcpy.AddMessage("\nOutput Maximum Upslope Grid: "+ssa)

drp=arcpy.GetParameterAsText(15)
if arcpy.Exists(drp):    
    arcpy.AddMessage("\nOutput Drop Analysis Table: "+drp)

usedroprange = arcpy.GetParameterAsText(16)
arcpy.AddMessage("\nSelect Threshold by Drop Analysis: "+usedroprange)

minthresh=arcpy.GetParameterAsText(17)
arcpy.AddMessage("\nMinimum Threshold Value: "+minthresh)

maxthresh=arcpy.GetParameterAsText(18)
arcpy.AddMessage("\nMaximum Threshold Value: "+maxthresh)

numthresh=arcpy.GetParameterAsText(19)
arcpy.AddMessage("\nNumber of Threshold Values: "+numthresh)

logspace=arcpy.GetParameterAsText(20)
arcpy.AddMessage("\nLogarithmic Spacing: "+logspace)

# Construct first command
cmd = 'mpiexec -n ' + inputProc + ' SlopeArea -slp ' + '"' + slp + '"' + ' -sca ' + '"' + sca + '"' + ' -sa ' + '"' + sa + '"' + ' -par ' + slpexp + ' ' + areaexp
arcpy.AddMessage("\nCommand Line: "+cmd)
# Submit command to operating system
os.system(cmd)
# Capture the contents of shell command and print it to the arcgis dialog box
process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
arcpy.AddMessage('\nProcess started:\n')
for line in process.stdout.readlines():
    arcpy.AddMessage(line)
arcpy.CalculateStatistics_management(sa)

# Construct second command
cmd = 'mpiexec -n ' + inputProc + ' D8FlowPathExtremeUp -p ' + '"' + p + '"' + ' -sa ' + '"' + sa + '"' + ' -ssa ' + '"' + ssa + '"'
if arcpy.Exists(shapelyr):
    cmd = cmd + ' -o ' + '"' + shapefile + '"'
##if maximumupslope == 'false':
##    cmd = cmd + ' -min '
if contcheck == 'false':
    cmd = cmd + ' -nc '

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
    cmd = cmd + ' -thresh ' + thresh    
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
