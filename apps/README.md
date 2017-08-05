Apps
====

This directory contains a number of useful apps which wrap subsets of RichDEM's
functionality for easy use.

The apps are as follows:

**rd_compare**: Determines whether, and in what ways, two rasters differ.
                Geotransform, NoData, Projection, Width, Height, and all data 
                values are checked.

**rd_geotransform**: Display or set the geotransform of a raster.

**rd_depressions_has**: Determines whether a raster contains depression(s).

**rd_depressions_mask**: Return a mask of all the depressions in a raster.
                         1 indicates cells that were in a depression,
                         0 indicates cells that were not in a depression,
                         3 indicates NoData.

**rd_depressions_flood**: Eliminate depressions by flooding them.

**rd_projection**: Alter a given raster's projection either by specifying it or
                   copying it from a second raster.

**rd_processing_history**: Display the processing history of a raster.

**rd_no_data**: Get or set a raster's NoData value

**rd_raster_display**: Prints a raster to the terminal

**rd_raster_inspect**: Print a subregion of a raster to the terminal

**rd_flow_accumulation**: Calculate flow accumulation in terms of upstream area
                          using one of a large number of algorithms.

**rd_terrain_property**: Calculate terrain properties such as slope, aspect, and
                         curvature.

TODO
====

rd_d8_flowdirs
rd_expand_dimensions
rd_flood_for_flowdirs
rd_hist
rd_layout_check.py
rd_layout_display.py
rd_layout_find_square.py
rd_loop_check
rd_merge_rasters_by_layout
rd_raster_to_tikz.py
rd_taudem_d8_to_richdem_d8
