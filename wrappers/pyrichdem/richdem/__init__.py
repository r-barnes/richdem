import pkg_resources
import _richdem
import datetime

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

def _AddAnalysis(arr, analysis):
  print("Add analysis: {0}".format(analysis))
  metastr  = "\n{nowdate} | {progname} v{progver} (hash={hash}, compiled={compdate}) | {analysis}".format(
    nowdate  = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S.%f UTC"),
    progname = "RichDEM (Python)",
    progver  = pkg_resources.require("richdem")[0].version,
    hash     = "Unknown",
    compdate = "Unknown",
    analysis = analysis
  )

  if not "PROCESSING_HISTORY" in arr.metadata:
    arr.metadata["PROCESSING_HISTORY"] = ""
  arr.metadata["PROCESSING_HISTORY"] += metastr


#2017-11-28 16:41:58 UTC | RichDEM v0.0.0 (hash=06ebc928f7542e84, compiled=2017-11-23 16:30:45 UTC) | rd_depressions_flood.exe /home/rick/data/gis/beauford.tif /z/out.tif 0

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
  ret.geotransform = _richdem.VecDouble(src_ds.GetGeoTransform())
  ret.setNoData(srcband.GetNoDataValue())
  for k,v in src_ds.GetMetadata().items():
    ret.metadata[k] = v

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
  for k,v in dem.metadata.items():
    data_set.SetMetadataItem(k,v)

def WrapNumPy(nparray):
  richdem_arrs = {
    'int8':    _richdem.Array2D_int8_t,
    'int16':   _richdem.Array2D_int16_t,
    'int32':   _richdem.Array2D_int32_t,
    'int64':   _richdem.Array2D_int64_t,
    'uint8':   _richdem.Array2D_uint8_t,
    'uint16':  _richdem.Array2D_uint16_t,
    'uint32':  _richdem.Array2D_uint32_t,
    'uint64':  _richdem.Array2D_uint64_t,
    'float32': _richdem.Array2D_float,
    'float64': _richdem.Array2D_double
  }

  if not nparray.dtype in richdem_arrs:
    raise Exception("No equivalent RichDEM datatype.")

  return richdem_arrs[nparray.dtype](nparray.data,nparray.shape[1],np.shape[0])

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
  _AddAnalysis(dem, "FillDepressions(dem, epsilon={0})".format(epsilon))
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

  if not method in facc_methods:
    raise Exception("Invalid FlowAccumulation method. Valid methods are: " + ', '.join(facc_methods.keys()))

  accum = _richdem.Array2D_double(dem, 0)

  _AddAnalysis(accum, "FlowAccumulation(dem, method={0})".format(method))

  facc_methods[method](dem,accum)

  return accum

def TerrainAttributes(
  dem,
  method,
  zscale = 1
):
  terrain_methods = {
    #"spi":                _richdem.TA_SPI,
    #"cti":                _richdem.TA_CTI,
    "slope_riserun":      _richdem.TA_slope_riserun,
    "slope_percentage":   _richdem.TA_slope_percentage,
    "slope_degrees":      _richdem.TA_slope_degrees,
    "slope_radians":      _richdem.TA_slope_radians,
    "aspect":             _richdem.TA_aspect,
    "curvature":          _richdem.TA_curvature,
    "planform_curvature": _richdem.TA_planform_curvature,
    "profile_curvature":  _richdem.TA_profile_curvature,
  }

  if not method in terrain_methods:
    raise Exception("Invalid TerrainAttributes method. Valid methods are: " + ', '.join(terrain_methods.keys()))

  result = _richdem.Array2D_float(dem, 0)

  _AddAnalysis(result, "TerrainAttributes(dem, method={0}, zscale={1})".format(method,zscale))

  terrain_methods[method](dem,result,zscale)

  return result