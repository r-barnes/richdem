# Script Name: RichDEM_ArcScript
# 
# Created By:  Richard Barnes
# Date:        5/10/2012

import arcpy
import os
import sys
import subprocess

PROGPATH="richdem"

def expand_layer_info(layerobj):
	arcpy.AddMessage("\n")
	arcpy.AddMessage("Base name: " + layerobj.baseName)
	arcpy.AddMessage("Path: " + layerobj.catalogPath)
#	arcpy.AddMessage("Children: " + layerobj.children)
#	arcpy.AddMessage("Children Expanded: " + layerobj.childrenExpanded)
	arcpy.AddMessage("Data element type: " + layerobj.dataElementType)
	arcpy.AddMessage("Data type: " + layerobj.dataType)
	arcpy.AddMessage("Extension: " + layerobj.extension)
	arcpy.AddMessage("File name: " + layerobj.file)
#	arcpy.AddMessage("Have full properties been retrieved?: " + layerobj.fullPropsRetrieved)
#	arcpy.AddMessage("Has the metadata been retrieved?: " + layerobj.metadataRetrieved)
	arcpy.AddMessage("User-assigned element name: " + layerobj.name)
	arcpy.AddMessage("File path: " + layerobj.path)


def main():
	inputDEM = arcpy.GetParameterAsText(0)
	desc = arcpy.Describe(inputDEM)

	expand_layer_info(desc)
	#inZfile=str(desc.catalogPath)
	arcpy.AddMessage("\nInput Elevation file: "+inputDEM)

	convert_dem_to_raster=arcpy.GetParameterAsText(1)
	arcpy.AddMessage("Convert DEM to Raster? " + convert_dem_to_raster)

	#import arcpy
	#from arcpy import env
	#env.workspace = "c:/data"
	#arcpy.RasterToASCII_conversion("elevation", "c:/output/sa500.asc")

	fill_depressions=arcpy.GetParameterAsText(2)
	arcpy.AddMessage("Fill depressions? " + fill_depressions)

	output_depression_filled_dem=arcpy.GetParameterAsText(3)
	arcpy.AddMessage("Output depression filled DEM? " + output_depression_filled_dem)

	use_dinf=arcpy.GetParameterAsText(4)
	arcpy.AddMessage("Use D-infinite flow metric? " + use_dinf)

	output_flowdirs_before_flat_resolution=arcpy.GetParameterAsText(5)
	arcpy.AddMessage("Output flowdirs before flat resolution? " + output_flowdirs_before_flat_resolution)

	output_flowdirs_after_flat_resolution=arcpy.GetParameterAsText(6)
	arcpy.AddMessage("Output flowdirs after flat resolution? " + output_flowdirs_after_flat_resolution)

	output_flow_acculm=arcpy.GetParameterAsText(7)
	arcpy.AddMessage("Output flow accumulation? " + output_flow_acculm)

	# Construct the taudem command line.  Put quotes around file names in case there are spaces

	#Z:\home\rick\projects\watershed\richdem\richdem  [-a <file>] [-f <file>]
	#                                        [-u <file>] [-l <file>] [-p] [-8]
	#                                        [--] [--version] [-h] <Input DEM>

	cmd = PROGPATH

	if fill_depressions=='true':
		cmd+=' -p'
		if output_depression_filled_dem:
			cmd+=' -l ' + output_depression_filled_dem

	if use_dinf!='true':
		cmd+=' -8 '

	if output_flowdirs_before_flat_resolution:
		cmd+=' -u ' + output_flowdirs_before_flat_resolution

	if output_flowdirs_after_flat_resolution:
		cmd+=' -f ' + output_flowdirs_after_flat_resolution

	if output_flow_acculm:
		cmd+=' -a ' + output_flow_acculm

	arcpy.AddMessage("\nCommand Line: "+cmd)
	os.system(cmd)
	process=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	arcpy.AddMessage('\nProcess started:\n')
	for line in process.stdout.readlines():
	    arcpy.AddMessage(line)

	#  Calculate statistics so that grids display with correct bounds
	#arcpy.AddMessage('Executing: Calculate Statistics\n')
	#arcpy.CalculateStatistics_management(outFile)

main()
