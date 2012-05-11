# Script Name: RichDEM_ArcScript
# 
# Created By:  Richard Barnes
# Date:        5/10/2012

import arcpy
import sys
import subprocess
import os
import tempfile
import uuid

PROGPATH="richdem.exe"

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

def make_temp_file():
	#temp_file_object=tempfile.NamedTemporaryFile(delete=False, suffix=".asc")
	#temp_file_name=temp_file_object.name
	#temp_file_object.close()
	temp_file_name=os.path.join(tempfile.gettempdir(),"rd"+str(uuid.uuid1())+".asc")
	arcpy.AddMessage("Created temporary file: " + temp_file_name)
	return temp_file_name

def main():
	global PROGPATH
	inputDEM = arcpy.GetParameterAsText(0)
	desc = arcpy.Describe(inputDEM)

	expand_layer_info(desc)
	arcpy.env.workspace=desc.path
	#inZfile=str(desc.catalogPath)
	arcpy.AddMessage("\nInput Elevation file: "+inputDEM)

	inputDEM_asc=make_temp_file()
	arcpy.RasterToASCII_conversion(inputDEM, inputDEM_asc)

	fill_depressions=arcpy.GetParameterAsText(1)
	arcpy.AddMessage("Fill depressions? " + fill_depressions)

	output_depression_filled_dem=arcpy.GetParameterAsText(2)
	arcpy.AddMessage("Output depression filled DEM? " + output_depression_filled_dem)

	use_dinf=arcpy.GetParameterAsText(3)
	arcpy.AddMessage("Use D-infinite flow metric? " + use_dinf)

	output_flow_acculm=arcpy.GetParameterAsText(4)
	arcpy.AddMessage("Output flow accumulation? " + output_flow_acculm)



	PROGPATH=os.path.join(sys.path[0],PROGPATH)
	cmd = [PROGPATH]

	output_depression_filled_dem_temp=""
	if fill_depressions=='true':
		cmd.append('-p')
		if output_depression_filled_dem:
			cmd.append('-l')
			output_depression_filled_dem_temp=make_temp_file()
			cmd.append(output_depression_filled_dem_temp)

	if use_dinf!='true':
		cmd.append('-8')

	output_flow_acculm_temp=""
	if output_flow_acculm:
		cmd.append('-a')
		output_flow_acculm_temp=make_temp_file()
		cmd.append(output_flow_acculm_temp)

	cmd.append(inputDEM_asc)

	arcpy.AddMessage("Command Line: "+str(cmd))
	try:
		process=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
	except:
		arcpy.AddMessage("Failed to run process!")
		sys.exit()
	arcpy.AddMessage('\nProcess started:\n')
	
	while True:
		line=process.stdout.readline()
		process.poll()
		if process.returncode!=None: break
		arcpy.AddMessage(line)

	process.poll()
	arcpy.AddMessage("\nProcess return code: " + str(process.returncode))
	if process.returncode!=0:
		arcpy.AddMessage("Process failed to run!")
		sys.exit()

	if fill_depressions=='true' and output_depression_filled_dem:
		arcpy.AddMessage('Converting Depression Filled DEM to Raster')
		arcpy.ASCIIToRaster_conversion(output_depression_filled_dem_temp, output_depression_filled_dem, "FLOAT")
		os.remove(output_depression_filled_dem_temp)

	if output_flow_acculm:
		arcpy.AddMessage('Converting Flow Accumulation DEM to Raster')
		arcpy.ASCIIToRaster_conversion(output_flow_acculm_temp, output_flow_acculm, "FLOAT")
		os.remove(output_flow_acculm_temp)

main()

#Note:
#A useful function is "arcpy.CalculateStatistics_management(FILE_NAME)"
