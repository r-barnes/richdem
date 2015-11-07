#!/bin/bash -x
mkdir /z/retain

#Generate a reference
time mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/data/beauford.tif /z/singlecore-

#Does block size make a difference?
time mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/data/beauford.tif /z/a500block-   -w 500 -h 500 
gdal_merge.py -o /z/500block.tif -of GTiff -n -9999 -a_nodata -9999 /z/a*tif

time mpirun -n 4 ./parallel_pit_fill.exe one @offloadall         ~/data/beauford.tif /z/b600block-   -w 600 -h 600 
gdal_merge.py -o /z/600block.tif -of GTiff -n -9999 -a_nodata -9999 /z/b*tif


#Does offloading make a difference?
time mpirun -n 50 ./parallel_pit_fill.exe one @retainall          ~/data/beauford.tif /z/c500block-retain-  -w 500 -h 500 
gdal_merge.py -o /z/500block-retain.tif -of GTiff -n -9999 -a_nodata -9999 /z/c500block-retain-beauford*tif

time mpirun -n 4  ./parallel_pit_fill.exe one /z/retain/500block- ~/data/beauford.tif /z/d500block-save-    -w 500 -h 500 
gdal_merge.py -o /z/500block-save.tif   -of GTiff -n -9999 -a_nodata -9999 /z/d500block-save-beauford*tif

#Calculate differences
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block.tif        --calc="A-B" --outfile=/z/500block-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/600block.tif        --calc="A-B" --outfile=/z/600block-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block-retain.tif --calc="A-B" --outfile=/z/500block-retain-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/singlecore-beauford-0-fill.tif -B /z/500block-save.tif   --calc="A-B" --outfile=/z/500block-save-diff.tif

gdalinfo -mm /z/500block-diff.tif        | grep Max
gdalinfo -mm /z/600block-diff.tif        | grep Max
gdalinfo -mm /z/500block-retain-diff.tif | grep Max
gdalinfo -mm /z/500block-save-diff.tif   | grep Max


#exit

cd /z/

echo "N41W084.hgt N41W085.hgt N41W086.hgt N42W084.hgt N42W085.hgt N42W086.hgt N43W084.hgt N43W085.hgt N43W086.hgt" | sed 's/ /\n/g' | awk '{print $1".zip"}' | xargs -n 1 -I {} unzip ~/temp/dds.cr.usgs.gov/srtm/version2_1/SRTM1/Region_03/{}

echo "N41W084.hgt,N41W085.hgt,N41W086.hgt
N42W084.hgt,N42W085.hgt,N42W086.hgt
N43W084.hgt,N43W085.hgt,N43W086.hgt" > layout

gdal_merge.py -o srtm_merged.tif -of GTiff -ot Int16 -n -32768 -a_nodata -32768 N41W084.hgt N41W085.hgt N41W086.hgt N42W084.hgt N42W085.hgt N42W086.hgt N43W084.hgt N43W085.hgt N43W086.hgt

cd -

time mpirun -n 2 ./parallel_pit_fill.exe one @offloadall   /z/srtm_merged.tif /z/srtm_merged_fill                                   #1m0.626s
time mpirun -n 4 ./parallel_pit_fill.exe one /z/retain/srtm_merged-   /z/srtm_merged.tif /z/srtm1000block-retain- -w 1000 -h 1000   #0m19.752s
time mpirun -n 4 ./parallel_pit_fill.exe one @offloadall              /z/srtm_merged.tif /z/srtm1000block-offload- -w 1000 -h 1000  #0m31.814s
time mpirun -n 4 ./parallel_pit_fill.exe many @offloadall             /z/layout          /z/srtm-sep-offload- -w 1000 -h 1000 -V -H #0m49.974s #TODO
time mpirun -n 4 ./parallel_pit_fill.exe many /z/retain/srtm_sep-     /z/layout          /z/srtm-sep-retain-  -w 1000 -h 1000 -V -H #0m19.462s

gdal_merge.py -o /z/output/srtm1000block-retain.tif  -of GTiff -n -9999 -a_nodata -9999 /z/srtm1000block-retain-*tif
gdal_merge.py -o /z/output/srtm1000block-offload.tif -of GTiff -n -9999 -a_nodata -9999 /z/srtm1000block-offload-*tif
gdal_merge.py -o /z/output/srtm-sep-offload.tif      -of GTiff -n -9999 -a_nodata -9999 /z/srtm-sep-offload-*tif
gdal_merge.py -o /z/output/srtm-sep-retain.tif       -of GTiff -n -9999 -a_nodata -9999 /z/srtm-sep-retain-*tif

gdal_calc.py --NoDataValue=-9999 -A /z/srtm_merged_fillsrtm_merged-0-fill.tif -B /z/output/srtm1000block-retain.tif  --calc="A-B" --outfile=/z/srtm1000block-retain-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/srtm_merged_fillsrtm_merged-0-fill.tif -B /z/output/srtm1000block-offload.tif --calc="A-B" --outfile=/z/srtm1000block-offload-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/srtm_merged_fillsrtm_merged-0-fill.tif -B /z/output/srtm-sep-offload.tif      --calc="A-B" --outfile=/z/srtm-sep-offload-diff.tif
gdal_calc.py --NoDataValue=-9999 -A /z/srtm_merged_fillsrtm_merged-0-fill.tif -B /z/output/srtm-sep-retain.tif       --calc="A-B" --outfile=/z/srtm-sep-retain-diff.tif

gdalinfo -mm /z/srtm1000block-retain-diff.tif  | grep Max
gdalinfo -mm /z/srtm1000block-offload-diff.tif | grep Max
gdalinfo -mm /z/srtm-sep-offload-diff.tif      | grep Max
gdalinfo -mm /z/srtm-sep-retain-diff.tif       | grep Max

#qgis 500block-retain.tif 500block-save.tif 500block.tif 600block.tif singlecore-beauford-0-fill.tif &