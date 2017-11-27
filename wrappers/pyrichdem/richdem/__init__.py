import _richdem
from _richdem import Array2D_float
from _richdem import Array2D_double
from _richdem import Array2D_uint8_t
from _richdem import Array2D_uint16_t
from _richdem import Array2D_uint32_t
from _richdem import Array2D_int8_t
from _richdem import Array2D_int16_t
from _richdem import Array2D_int32_t

try:
  import numpy as np
  NUMPY_AVAILABLE = True
except:
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

  gdal_to_richdem = {
    gdal.GDT_Byte:    _richdem.Array2D_uint8_t,
    gdal.GDT_Int16:   _richdem.Array2D_int16_t,
    gdal.GDT_Int32:   _richdem.Array2D_int32_t,
    gdal.GDT_UInt16:  _richdem.Array2D_uint16_t,
    gdal.GDT_UInt32:  _richdem.Array2D_uint32_t,
    gdal.GDT_Float32: _richdem.Array2D_float,
    gdal.GDT_Float64: _richdem.Array2D_double
  }


  #Read in data
  src_ds     = gdal.Open(filename)
  srcband    = src_ds.GetRasterBand(1)
  srcdata    = srcband.ReadAsArray()
  # raster_srs = osr.SpatialReference()
  # raster_srs.ImportFromWkt(raster.GetProjectionRef())

  if not srcband.DataType in gdal_to_richdem:
    raise Exception("This datatype is not supported. Please file a bug report on RichDEM.")

  ret = gdal_to_richdem[srcband.DataType]()
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


def FillDepressions(
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

def FlowAccumulation(
  dem,
  method = 'D8'
):
  facc_methods = {
    "Tarboton":          _richdem.FA_Tarboton,
    "Dinf":              _richdem.FA_Tarboton,
    "Holmgren":          _richdem.FA_Holmgren,
    "Quinn":             _richdem.FA_Quinn,
    "Freeman":           _richdem.FA_Freeman,
    "FairfieldLeymarie": _richdem.FA_FairfieldLeymarie,
    "Rho8":              _richdem.FA_Rho8,
    "OCallaghan":        _richdem.FA_OCallaghan,
    "D8":                _richdem.FA_D8
  }

  accum = _richdem.Array2D_double()

  facc_methods[method](dem,accum)

  return accum