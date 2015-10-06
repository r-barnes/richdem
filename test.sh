#!/bin/bash
mkdir /z/retain

#Generate a reference
mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/projects/watershed/data/beauford.tif /z/singlecore-

#Does block size make a difference?
mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/projects/watershed/data/beauford.tif /z/a500block-   -w 500 -h 500 
mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/projects/watershed/data/beauford.tif /z/b600block-   -w 600 -h 600 
gdal_merge.py -o /z/500block.tif -of GTiff -n -9999 -a_nodata -9999 /z/a*tif
gdal_merge.py -o /z/600block.tif -of GTiff -n -9999 -a_nodata -9999 /z/b*tif


#Does offloading make a difference?
mpirun -n 50 ./parallel_pit_fill.exe one @retainall          ~/projects/watershed/data/beauford.tif /z/c500block-retain-  -w 500 -h 500 
mpirun -n 4  ./parallel_pit_fill.exe one /z/retain/500block- ~/projects/watershed/data/beauford.tif /z/d500block-save-    -w 500 -h 500 

gdal_merge.py -o /z/500block-retain.tif -of GTiff -n -9999 -a_nodata -9999 /z/c500block-retain-beauford*tif
gdal_merge.py -o /z/500block-save.tif   -of GTiff -n -9999 -a_nodata -9999 /z/d500block-save-beauford*tif

gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block.tif        --calc="A-B" --outfile=/z/500block-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/600block.tif        --calc="A-B" --outfile=/z/600block-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block-retain.tif --calc="A-B" --outfile=/z/500block-retain-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block-save.tif   --calc="A-B" --outfile=/z/500block-save-diff.tif

gdalinfo -mm /z/500block-diff.tif        | grep Max
gdalinfo -mm /z/600block-diff.tif        | grep Max
gdalinfo -mm /z/500block-retain-diff.tif | grep Max
gdalinfo -mm /z/500block-save-diff.tif   | grep Max

#qgis 500block-retain.tif 500block-save.tif 500block.tif 600block.tif singlecore-beauford-0-fill.tif &