import _richdem
from _richdem import Array2Dfloat

try:
  import numpy as np
  NUMPY_AVAILABLE = True
else:
  NUMPY_AVAILABLE = False

try:
  from osgeo import gdal
  gdal.UseExceptions()
  GDAL_AVAILABLE = True
except:
  GDAL_AVAILABLE = False

def LoadGDAL(filename):
  if not GDAL_AVAILABLE:
    raise Exception("richdem.LoadGDAL() requires GDAL.")

  #Read in data
  src_ds     = gdal.Open(filename)
  srcband    = src_ds.GetRasterBand(1)
  srcdata    = srcband.ReadAsArray()
  # raster_srs = osr.SpatialReference()
  # raster_srs.ImportFromWkt(raster.GetProjectionRef())

  ret = Array2Dfloat()
  ret.fromArray(srcdata)
  ret.projection   = src_ds.GetProjectionRef()
  ret.geotransform = src_ds.GetGeoTransform()
  ret.setNoData(srcband.GetNoDataValue())

  print(src_ds.GetMetadata())

  return ret

def SaveGDAL(filename, dem):
  if not GDAL_AVAILABLE:
    raise Exception("richdem.SaveGDAL() requires GDAL.")
  if not NUMPY_AVAILABLE:
    raise Exception("richdem.SaveGDAL() requires NumPy.")

  driver    = gdal.GetDriverByName('GTiff')
  data_type = gdal.GDT_Float32 #TODO
  data_set  = driver.Create(filename, xsize=dem.width(), ysize=dem.height(), bands=1, eType=data_type)
  data_set.SetGeoTransform(dem.geotransform)
  data_set.SetProjection(dem.projection)
  band = data_set.GetRasterBand(1)
  band.SetNoDataValue(dem.noData())
  band.WriteArray(np.array(dem))


def fillDepressions(
  dem,
  epsilon = False,
):
  """Fills all depressions in an elevation model.

     Parameters:
     dem -- A NumPy elevation model

     Returns:
     Modified dem in-place to remove depressions
  """
  if epsilon:
    return _richdem.rdPFepsilon(dem)
  else:
    return _richdem.rdFillDepressions(dem)
