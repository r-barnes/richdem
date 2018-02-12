Title: Efficient Priority-Flood depression filling in raster digital elevation models

This repository contains the source codes of the algorithm presented in the manuscript above. These codes were used in performing the tests described in the manuscript.

The coded can be compiled on both Windows and Linux platforms.

FillDEM supports floating-point GeoTIFF file format through the GDAL library. Please include GDAL library into your compilation. GDAL 1.91 was used in our experiments.

The algorithms described in the manuscript can be run using the following command line:

fill fillDEM Input_DEM Output_DEM

Example: ./fill fillDEM dem.tif. dem_filled.tif

Use the algorithm proposed in the manuscript to fill the input DEM "dem.tif" and create a filled DEM "dem_filled.tif".

For check the correctness of our algorithm easily, this program provide a function that can be used to generate synthetic DEM.

This function can be run using the following command line:

fill perlinDEM rows cols dem.tif

Example: ./fill perlinDEM 100 100 dem.tif
 
The test data used in the manuscript can be downloaded at http://www.mngeo.state.mn.us. You need ArcGIS to convert the DEM into GeoTIFF format.
