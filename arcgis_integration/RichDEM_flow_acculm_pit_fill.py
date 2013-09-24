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
import time

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
	#temp_file_object=tempfile.NamedTemporaryFile(delete=False)
	#temp_file_name=temp_file_object.name
	#temp_file_object.close()
	temp_file_name=os.path.join(tempfile.gettempdir(),"rd"+str(uuid.uuid1()))
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

	start = time.time()
	inputDEM_flt=make_temp_file()
	arcpy.RasterToFloat_conversion(inputDEM, inputDEM_flt+".flt")
	arcpy.AddMessage("Creating the temporary file took and raster conversion took " + str(round(time.time() - start,2)) + "s.")

	fill_depressions=arcpy.GetParameterAsText(1)
	arcpy.AddMessage("Fill depressions? " + fill_depressions)

	output_depression_filled_dem=arcpy.GetParameterAsText(2)
	arcpy.AddMessage("Output depression filled DEM? " + output_depression_filled_dem)

	use_dinf=arcpy.GetParameterAsText(3)
	arcpy.AddMessage("Use D-infinite flow metric? " + use_dinf)

	output_flow_acculm=arcpy.GetParameterAsText(4)
	arcpy.AddMessage("Output flow accumulation? " + output_flow_acculm)
	
	zscale=arcpy.GetParameterAsText(5)
	arcpy.AddMessage("z-scaling factor: " + zscale)
	
	output_cti=arcpy.GetParameterAsText(6)
	arcpy.AddMessage("Output CTI? " + output_cti)
	
	output_spi=arcpy.GetParameterAsText(7)
	arcpy.AddMessage("Output SPI? " + output_spi)



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
		
	if zscale:
		cmd.append('--zscale')
		cmd.append(zscale)
		
	if output_cti:
		cmd.append('--cti')
		output_cti_temp=make_temp_file()
		cmd.append(output_cti_temp)
		
	if output_spi:
		cmd.append('--spi')
		output_spi_temp=make_temp_file()
		cmd.append(output_spi_temp)

	cmd.append(inputDEM_flt+".flt")

	arcpy.AddMessage("Command Line: "+str(cmd))
	try:
		startupinfo=subprocess.STARTUPINFO()
		startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
		process=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, startupinfo=startupinfo)
	except:
		arcpy.AddMessage("Failed to run process!")
		sys.exit()
	arcpy.AddMessage('\nProcess started:\n')

	while True:
		process.poll()
		if process.returncode!=None:
			break
		line=process.stdout.readline().rstrip()
		if line=='':
			continue

		if line[0]=='%':
			arcpy.ResetProgressor()
			arcpy.SetProgressor("step", line[1:])
			arcpy.AddMessage(line[1:])
		elif line[0:3]=='P%c':
			arcpy.SetProgressorPosition(0)
		elif line[0:2]=='P%':
			arcpy.SetProgressorPosition(int(line[2:]))
		else:
			arcpy.AddMessage(line)

	process.poll()
	arcpy.AddMessage("\nProcess return code: " + str(process.returncode))
	if process.returncode!=0:
		arcpy.AddMessage("Process failed to run!")
		sys.exit()

	if fill_depressions=='true' and output_depression_filled_dem:
		arcpy.AddMessage('Converting Depression Filled to Raster')
		arcpy.FloatToRaster_conversion(output_depression_filled_dem_temp+".flt", output_depression_filled_dem)
		os.remove(output_depression_filled_dem_temp+".flt")

	if output_flow_acculm:
		arcpy.AddMessage('Converting Flow Accumulation to Raster')
		arcpy.FloatToRaster_conversion(output_flow_acculm_temp+".flt", output_flow_acculm)
		os.remove(output_flow_acculm_temp+".flt")
		
	if output_cti:
		arcpy.AddMessage('Converting CTI to Raster')
		arcpy.FloatToRaster_conversion(output_cti_temp+".flt", output_cti)
		os.remove(output_cti_temp+".flt")
		
	if output_spi:
		arcpy.AddMessage('Converting SPI to Raster')
		arcpy.FloatToRaster_conversion(output_spi_temp+".flt", output_spi)
		os.remove(output_spi_temp+".flt")

main()

#Note:
#A useful function is "arcpy.CalculateStatistics_management(FILE_NAME)"
