# Script Name: WatershedGridToShapefile
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
w=str(desc.catalogPath)
arcpy.AddMessage("\nInput Elevation file: "+w)

# Output
shfl = arcpy.GetParameterAsText(1)
arcpy.AddMessage("\nOutput Stream Source file: "+shfl)

# Convert tiff to shp
cmd = arcpy.RasterToPolygon_conversion(w, shfl, "NO_SIMPLIFY", "Value")
