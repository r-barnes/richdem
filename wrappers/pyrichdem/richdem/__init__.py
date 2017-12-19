import pkg_resources
import _richdem
import datetime
import copy

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

def _RichDEMVersion():
  return "RichDEM (Python {pyver}) (hash={hash}, hashdate={compdate})".format(
    pyver    = pkg_resources.require("richdem")[0].version,
    hash     = _richdem.rdHash(),
    compdate = _richdem.rdCompileTime()
  )

def _AddAnalysis(arr, analysis):
  metastr  = "\n{nowdate} | {verstr} | {analysis}".format(
    nowdate  = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S.%f UTC"),
    verstr   = _RichDEMVersion(),
    analysis = analysis
  )

  if not "PROCESSING_HISTORY" in arr.metadata:
    arr.metadata["PROCESSING_HISTORY"] = ""
  arr.metadata["PROCESSING_HISTORY"] += metastr

def rdShow(rda, ignore_colours=[], show=True):
  try:
    import matplotlib.pyplot as plt
    import matplotlib
  except:
    raise Exception("matplotlib must be installed to use rdShow!")
  disparr = np.array(rda, copy=True)
  disparr[disparr==rda.no_data] = np.nan
  for c in ignore_colours:
    disparr[disparr==c] = np.nan
  vmin, vmax = np.nanpercentile(disparr, [2, 98])
  #current_cmap = matplotlib.cm.get_cmap()
  #current_cmap.set_bad(color='red')
  plt.imshow(disparr, vmin=vmin, vmax=vmax)
  plt.set_cmap('gray')
  if show:
    plt.show()


class rdarray(np.ndarray):
  def __new__(cls, array, meta_obj=None, dtype=None, order=None, **kwargs):
    obj = np.asarray(array, dtype=dtype, order=order).view(cls) 
    
    if meta_obj is not None:
      obj.metadata     = copy.deepcopy(getattr(meta_obj, 'metadata',     None))
      obj.no_data      = copy.deepcopy(getattr(meta_obj, 'no_data',      None))
      obj.projection   = copy.deepcopy(getattr(meta_obj, 'projection',   None))
      obj.geotransform = copy.deepcopy(getattr(meta_obj, 'geotransform', None))

    return obj

  def __array_finalize__(self, obj):
    if obj is None: return
    self.metadata     = copy.deepcopy(getattr(obj, 'metadata',     None))
    self.no_data      = copy.deepcopy(getattr(obj, 'no_data',      None))
    self.projection   = copy.deepcopy(getattr(obj, 'projection',   None))
    self.geotransform = copy.deepcopy(getattr(obj, 'geotransform', None))

  def wrap(self):
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
    dtype = str(self.dtype)
    if not dtype in richdem_arrs:
      raise Exception("No equivalent RichDEM datatype.")
    rda = richdem_arrs[dtype](self)
    rda.setNoData(self.no_data)
    if self.geotransform:
      rda.geotransform = np.array(self.geotransform, dtype='float64')
    else:
      print("Warning! No geotransform defined. Choosing a standard one!")
      rda.geotransform = np.array([0,1,0,0,0,-1], dtype='float64')
    return rda 

  def copyFromWrapped(self, wrapped):
    self.no_data   = wrapped.noData()
    #print("IS IT IN",("PROCESSING_HISTORY" in wrapped.metadata))
    #self.metadata += "\n"+wrapped.metadata["PROCESSING_HISTORY"].replace("\n","\t\n")


def LoadGDAL(filename, no_data=None):
  """Read a GDAL file.

     Opens any file GDAL can read, selects the first raster band, and loads it
     and its metadata into a RichDEM array of the appropriate data type.

     If you need to do something more complicated, look at the source of this
     function.

     Parameters:
     filename -- Name of the raster file to open
     no_data  -- Optionally, set the no_data value to this.

     Returns:
     A RichDEM array
  """  
  if not GDAL_AVAILABLE:
    raise Exception("richdem.LoadGDAL() requires GDAL.")

  allowed_types = {gdal.GDT_Byte,gdal.GDT_Int16,gdal.GDT_Int32,gdal.GDT_UInt16,gdal.GDT_UInt32,gdal.GDT_Float32,gdal.GDT_Float64}

  #Read in data
  src_ds  = gdal.Open(filename)
  srcband = src_ds.GetRasterBand(1)
  srcdata = rdarray(srcband.ReadAsArray())

  # raster_srs = osr.SpatialReference()
  # raster_srs.ImportFromWkt(raster.GetProjectionRef())

  if not srcband.DataType in allowed_types:
    raise Exception("This datatype is not supported. Please file a bug report on RichDEM.")

  srcdata.projection   = src_ds.GetProjectionRef()
  srcdata.geotransform = src_ds.GetGeoTransform()

  if no_data is None:
    srcband_no_data = srcband.GetNoDataValue()
    if srcband_no_data is None:
      raise Exception("The source data did not have a NoData value. Please use the no_data argument to specify one. If should not be equal to any of the actual data values. If you are using all possible data values, then the situation is pretty hopeless - sorry.")
    else:
      srcdata.no_data = srcband.GetNoDataValue()
  else:
    srcdata.no_data = no_data

  srcdata.metadata = dict()
  for k,v in src_ds.GetMetadata().items():
    srcdata.metadata[k] = v

  _AddAnalysis(srcdata, "LoadGDAL(filename={0}, no_data={1})".format(filename, no_data))

  return srcdata



def SaveGDAL(filename, dem):
  """Save a GDAL file.

     Saves a RichDEM array to a data file in GeoTIFF format.

     If you need to do something more complicated, look at the source of this
     function.

     Parameters:
     filename -- Name of the raster file to be created
     dem      -- RichDEM array to save.

     Returns:
     No Return
  """    
  if not GDAL_AVAILABLE:
    raise Exception("richdem.SaveGDAL() requires GDAL.")
  if not NUMPY_AVAILABLE:
    raise Exception("richdem.SaveGDAL() requires NumPy.")

  driver    = gdal.GetDriverByName('GTiff')
  data_type = gdal.GDT_Float32 #TODO
  data_set  = driver.Create(filename, xsize=dem.shape[1], ysize=dem.shape[0], bands=1, eType=data_type)
  data_set.SetGeoTransform(dem.geotransform)
  data_set.SetProjection(dem.projection)
  band = data_set.GetRasterBand(1)
  band.SetNoDataValue(dem.no_data)
  band.WriteArray(np.array(dem))
  for k,v in dem.metadata.items():
    data_set.SetMetadataItem(str(k),str(v))



def FillDepressions(
  dem,
  epsilon  = False,
  in_place = False
):
  """Fills all depressions in a DEM using in-place modification.

     Parameters:
     dem     -- An elevation model
     epsilon -- If True, an epsilon gradient is imposed to all flat regions.
                This ensures that there is always a local gradient.

     Returns:
     DEM without depressions.
  """

  if not in_place:
    print("Copying dem")
    dem = dem.copy()

  _AddAnalysis(dem, "FillDepressions(dem, epsilon={0})".format(epsilon))

  demw = dem.wrap()

  if epsilon:
    _richdem.rdPFepsilon(demw)
  else:
    _richdem.rdFillDepressions(demw)

  dem.copyFromWrapped(demw)

  if not in_place:
    print("Returning dem")
    return dem



def BreachDepressions(
  dem,
  mode           = "Complete",
  fill           = False,
  max_path_len   = None,
  max_path_depth = None,
  in_place       = False
):
  """Attempts to breach depressions in a DEM using in-place modification.

     Parameters:
     dem            -- An elevation model
     mode           -- String. Breaching mode to use (see below).
     fill           -- True/False. If depressions can't be breached, 
                       should they be filled?
     max_path_len   -- Maximum path length in cells that can be dug while
                       breaching a depression.
     max_path_depth -- Maximum depth in z-units that can be dug along the 
                       breaching path.

     Modes:
     Complete:    Breach everything.
                  Ignore max_path_len, max_path_depth.
                  There will be no depressions.
                  There will be no mercy.
     Selective:   Only breach those depressions that can be breached using the
                  above criteria.
     Constrained: Dig as long a path as necessary, but don't dig it deeper than
                  max_path_depth.

     Returns:
     DEM with depressions breached, filled, or left unaltered, per the above
     options.
  """

  if not in_place:
    dem = dem.copy()

  _AddAnalysis(dem, "BreachDepressions(dem, mode={0}, fill={1}, max_path_len={2}, max_path_depth={3})".format(mode, fill, max_path_len, max_path_depth))

  modes = {
    "Complete":    1,
    "Selective":   2,
    "Constrained": 3,
  }

  if not mode in modes:
    raise Exception("Unrecognised mode. Valid choices are: " + ', '.join(modes.keys()))

  if mode!="Complete" and max_path_len is None:
    raise Exception("Must provide a 'max_path_len'")
  elif mode!="Complete" and max_path_depth is None:
    raise Exception("Must provide a 'max_path_depth'")
  elif mode=="Complete":
    max_path_len   = 0
    max_path_depth = 0

  demw = dem.wrap()

  _richdem.rdBreach(demw, modes[mode], fill, max_path_len, max_path_depth)

  dem.copyFromWrapped(demw)

  if not in_place:
    return dem



def FlowAccumulation(
  dem,
  method   = None,
  exponent = None
):
  """Calculates flow accumulation. A variety of methods are available.

     Parameters:
     dem      -- An elevation model
     method   -- Flow accumulation method to use. (See below.)
     exponent -- Some methods require an exponent; refer to the relevant
                 publications for details.

     Method            Note                           Reference
     Tarboton          Alias for Dinf.
     Dinf              Alias for Tarboton.
     Quinn             Holmgren with exponent=1.
     Holmgren(E)       Generalization of Quinn.
     Freeman(E)        TODO
     FairfieldLeymarie Alias for Rho8.
     Rho8              Alias for FairfieldLeymarie.
     OCallaghan        Alias for D8.                  10.1016/S0734-189X(84)80011-0
     D8                Alias for OCallaghan.          10.1016/S0734-189X(84)80011-0

     Methods marked (E) require the exponent argument.

     Returns:
     Flow accumulation according to the desired method.
  """
  facc_methods = {
    "Tarboton":          _richdem.FA_Tarboton,
    "Dinf":              _richdem.FA_Tarboton,
    "Quinn":             _richdem.FA_Quinn,
    "FairfieldLeymarie": _richdem.FA_FairfieldLeymarie,
    "Rho8":              _richdem.FA_Rho8,
    "OCallaghan":        _richdem.FA_OCallaghan,
    "D8":                _richdem.FA_D8
  }

  facc_methods_exponent = {
    "Freeman":           _richdem.FA_Freeman,
    "Holmgren":          _richdem.FA_Holmgren
  }

  accum  = rdarray(np.zeros(shape=dem.shape, dtype='float64'), meta_obj=dem)
  accumw = accum.wrap()

  _AddAnalysis(accum, "FlowAccumulation(dem, method={0})".format(method))

  if method in facc_methods:
    facc_methods[method](dem.wrap(),accumw)
  elif method in facc_methods_exponent:
    if exponent is None:
      raise Exception('FlowAccumulation method "'+method+'" requires an exponent!')
    facc_methods_exponent[method](dem.wrap(),accumw,exponent)
  else:
    raise Exception("Invalid FlowAccumulation method. Valid methods are: " + ', '.join(list(facc_methods.keys()) + list(facc_methods_exponent.keys()) ))

  accum.copyFromWrapped(accumw)

  return accum


def TerrainAttribute(
  dem,
  attrib,
  zscale = 1.0
):
  """Calculates terrain attributes. A variety of methods are available.

     Parameters:
     dem      -- An elevation model
     attrib   -- Terrain attribute to calculate. (See below.)
     zscale   -- How much to scale the z-axis by prior to calculation

     Method:
     slope_riserun
     slope_percentage
     slope_degrees
     slope_radians
     aspect
     curvature
     planform_curvature
     profile_curvature

     Returns:
     A raster of the indicated terrain attributes.
  """
  terrain_attribs = {
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

  if not attrib in terrain_attribs:
    raise Exception("Invalid TerrainAttributes attribute. Valid attributes are: " + ', '.join(terrain_attribs.keys()))

  result  = rdarray(np.zeros(shape=dem.shape, dtype='float32'), meta_obj=dem)
  resultw = result.wrap()

  _AddAnalysis(result, "TerrainAttribute(dem, attrib={0}, zscale={1})".format(attrib,zscale))

  terrain_attribs[attrib](dem.wrap(),resultw,zscale)

  result.copyFromWrapped(resultw)

  return result
