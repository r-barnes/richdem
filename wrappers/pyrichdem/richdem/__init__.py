import copy
import datetime
import pkg_resources
from typing import Any, Dict, Final, List, Iterable, Optional, Tuple, Union

import numpy as np

try:
    import _richdem
except ImportError as e:
    print("COULD NOT LOAD RichDEM ENGINE! NOTHING WILL WORK!")
    raise e

from _richdem import depression_hierarchy

try:
    from osgeo import gdal

    gdal.UseExceptions()
    GDAL_AVAILABLE = True
except:
    GDAL_AVAILABLE = False

STANDARD_GEOTRANSFORM: Final[np.ndarray] = np.array([0, 1, 0, 0, 0, -1])

def _RichDEMVersion() -> str:
    return "RichDEM (Python {pyver}) (hash={hash}, hashdate={compdate})".format(
        pyver=pkg_resources.require("richdem")[0].version,
        hash=_richdem.rdHash(),
        compdate=_richdem.rdCompileTime(),
    )


def _AddAnalysis(rda: "rdarray", analysis: str) -> None:
    if type(rda) not in {rdarray, rd3array}:
        raise Exception("An rdarray or rd3array is required!")

    metastr = "\n{nowdate} | {verstr} | {analysis}".format(
        nowdate=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S.%f UTC"),
        verstr=_RichDEMVersion(),
        analysis=analysis,
    )

    if rda.metadata is None:
        rda.metadata = dict()
    if "PROCESSING_HISTORY" not in rda.metadata:
        rda.metadata["PROCESSING_HISTORY"] = ""
    rda.metadata["PROCESSING_HISTORY"] += metastr


def rdShow(
    rda: "rdarray",
    ignore_colours=[],
    show: bool = True,
    axes: bool = True,
    cmap: str = "gray",
    vmin: Optional[int] = None,
    vmax: Optional[int] = None,
    xmin: Optional[int] = None,
    xmax: Optional[int] = None,
    ymin: Optional[int] = None,
    ymax: Optional[int] = None,
    zxmin: Optional[int] = None,
    zxmax: Optional[int] = None,
    zymin: Optional[int] = None,
    zymax: Optional[int] = None,
    figsize: Tuple[int, int] = (4, 4),
    zcolor: str = "red",
    zloc: int = 1,
):
    if type(rda) is np.ndarray:
        rda = rdarray(rda)
    elif type(rda) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    try:
        import matplotlib.pyplot as plt
        import matplotlib
    except:
        raise Exception("matplotlib must be installed to use rdShow!")

    zoom_vars = [zxmin, zxmax, zymin, zymax]
    some_zoom = any(x is not None for x in zoom_vars)
    all_zoom = all(x is not None for x in zoom_vars)

    if some_zoom and not all_zoom:
        raise Exception("All zoom limits must be set for zooming to work!")
    elif all_zoom:
        try:
            # from mpl_toolkits.axes_grid1.inset_locator import mark_inset
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            from matplotlib.patches import Rectangle
        except:
            raise Exception("mpl_toolkits.axes_grid1 must be available!")

    disparr = np.array(rda, copy=True)
    disparr[disparr == rda.no_data] = np.nan
    for c in ignore_colours:
        disparr[disparr == c] = np.nan
    vmin_calc, vmax_calc = np.nanpercentile(disparr, [2, 98])
    if vmin is None:
        vmin = vmin_calc
    if vmax is None:
        vmax = vmax_calc

    fig, (ax, cax) = plt.subplots(
        ncols=2, figsize=figsize, gridspec_kw={"width_ratios": [1, 0.05]}
    )

    # current_cmap = matplotlib.cm.get_cmap()
    # current_cmap.set_bad(color='red')
    iax = ax.imshow(disparr, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    if all_zoom:
        assert zxmin is not None
        assert zxmax is not None
        assert zymin is not None
        assert zymax is not None
        axins = inset_axes(
            ax, width=2, height=2, loc=zloc, borderpad=0
        )  # , bbox_to_anchor=(0.9, -0.05, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        axins.set_xlim(xmin=zxmin, xmax=zxmax)
        axins.set_ylim(ymin=zymin, ymax=zymax)
        plt.setp(axins.get_xticklabels(), visible=False)
        plt.setp(axins.get_yticklabels(), visible=False)
        plt.setp(axins.get_xticklines(), visible=False)
        plt.setp(axins.get_yticklines(), visible=False)
        # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec='r', lw=2, ls=None) #ec='1.0' = top of colormap
        ax.add_patch(
            Rectangle(
                (zxmin, zymin),
                zxmax - zxmin,
                zymax - zymin,
                facecolor="none",
                edgecolor=zcolor,
                lw=2,
            )
        )  # , transform=fig.transFigure, facecolor='none'))
        plt.setp(tuple(axins.spines.values()), color=zcolor, lw=2)
        axins.imshow(disparr, vmin=vmin, vmax=vmax, cmap=cmap)

    fig.colorbar(iax, cax=cax)

    plt.tight_layout()

    if not axes:
        ax.axis("off")
    if show:
        plt.show()
    return {"vmin": vmin, "vmax": vmax}


class rdarray(np.ndarray):
    def __new__(
        cls, array, meta_obj=None, no_data: Optional[Union[float, int]]=None, dtype=None, order=None, geotransform: Optional[Iterable[float]]=None, **kwargs: Any
    ) -> "rdarray":
        obj = np.asarray(array, dtype=dtype, order=order).view(cls)

        if meta_obj is not None:
            obj.metadata = copy.deepcopy(getattr(meta_obj, "metadata", dict()))
            obj.no_data = copy.deepcopy(getattr(meta_obj, "no_data", None))
            obj.projection = copy.deepcopy(getattr(meta_obj, "projection", ""))
            obj.geotransform = copy.deepcopy(getattr(meta_obj, "geotransform", None))
        elif geotransform is not None:
            obj.geotransform = geotransform

        if no_data is not None:
            obj.no_data = no_data

        if no_data is None:
            raise Exception("A no_data value must be specified!")

        return obj

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return
        self.metadata = copy.deepcopy(getattr(obj, "metadata", dict()))
        self.no_data = copy.deepcopy(getattr(obj, "no_data", None))
        self.projection = copy.deepcopy(getattr(obj, "projection", ""))
        self.geotransform = copy.deepcopy(getattr(obj, "geotransform", None))

    def wrap(self):
        richdem_arrs = {
            "int8": _richdem.Array2D_int8_t,
            "int16": _richdem.Array2D_int16_t,
            "int32": _richdem.Array2D_int32_t,
            "int64": _richdem.Array2D_int64_t,
            "uint8": _richdem.Array2D_uint8_t,
            "uint16": _richdem.Array2D_uint16_t,
            "uint32": _richdem.Array2D_uint32_t,
            "uint64": _richdem.Array2D_uint64_t,
            "float32": _richdem.Array2D_float,
            "float64": _richdem.Array2D_double,
        }
        dtype = str(self.dtype)
        if dtype not in richdem_arrs:
            raise Exception(f"No equivalent RichDEM datatype to '{dtype}'.")

        rda = richdem_arrs[dtype](self)

        if self.no_data is None:
            print("Warning! no_data was None. Setting it to -9999!")
            rda.setNoData(-9999)
        else:
            rda.setNoData(self.no_data)

        if self.geotransform is not None:
            rda.geotransform = np.array(self.geotransform, dtype="float64")
        else:
            print(
                "Warning! No geotransform defined. Choosing a standard one! (Top left cell's top let corner at <0,0>; cells are 1x1.)"
            )
            rda.geotransform = np.array([0, 1, 0, 0, 0, -1], dtype="float64")

        return rda

    def copyFromWrapped(self, wrapped):
        self.no_data = wrapped.noData()
        # print("IS IT IN",("PROCESSING_HISTORY" in wrapped.metadata))
        # self.metadata += "\n"+wrapped.metadata["PROCESSING_HISTORY"].replace("\n","\t\n")


class rd3array(np.ndarray):
    def __new__(cls, array, meta_obj=None, no_data=None, order=None, **kwargs: Any) -> "rd3array":
        obj = np.asarray(array, dtype=np.float32, order=order).view(cls)

        if meta_obj is not None:
            obj.metadata = copy.deepcopy(getattr(meta_obj, "metadata", dict()))
            obj.no_data = copy.deepcopy(getattr(meta_obj, "no_data", None))
            obj.projection = copy.deepcopy(getattr(meta_obj, "projection", ""))
            obj.geotransform = copy.deepcopy(getattr(meta_obj, "geotransform", None))

        if no_data is not None:
            obj.no_data = no_data

        if no_data is None:
            raise Exception("A no_data value must be specified!")

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.metadata = copy.deepcopy(getattr(obj, "metadata", dict()))
        self.no_data = copy.deepcopy(getattr(obj, "no_data", None))
        self.projection = copy.deepcopy(getattr(obj, "projection", ""))
        self.geotransform = copy.deepcopy(getattr(obj, "geotransform", None))

    def wrap(self):
        richdem_arrs = {"float32": _richdem.Array3D_float}
        dtype = str(self.dtype)
        if dtype not in richdem_arrs:
            raise Exception("No equivalent RichDEM datatype.")

        rda = richdem_arrs[dtype](self)

        if self.no_data is None:
            print("Warning! no_data was None. Setting it to -9999!")
            rda.setNoData(-9999)
        else:
            rda.setNoData(self.no_data)

        if self.geotransform:
            rda.geotransform = np.array(self.geotransform, dtype="float64")
        else:
            print(
                "Warning! No geotransform defined. Choosing a standard one! (Top left cell's top let corner at <0,0>; cells are 1x1.)"
            )
            rda.geotransform = np.array([0, 1, 0, 0, 0, -1], dtype="float64")

        return rda

    def copyFromWrapped(self, wrapped):
        self.no_data = wrapped.noData()
        # print("IS IT IN",("PROCESSING_HISTORY" in wrapped.metadata))
        # self.metadata += "\n"+wrapped.metadata["PROCESSING_HISTORY"].replace("\n","\t\n")


def LoadGDAL(filename: str, no_data: Optional[float] = None) -> rdarray:
    """Read a GDAL file.

    Opens any file GDAL can read, selects the first raster band, and loads it
    and its metadata into a RichDEM array of the appropriate data type.

    If you need to do something more complicated, look at the source of this
    function.

    Args:
        filename (str):    Name of the raster file to open
        no_data  (float):  Optionally, set the no_data value to this.

    Returns:
        A RichDEM array
    """
    if not GDAL_AVAILABLE:
        raise Exception("richdem.LoadGDAL() requires GDAL.")

    allowed_types = {
        gdal.GDT_Byte,
        gdal.GDT_Int16,
        gdal.GDT_Int32,
        gdal.GDT_UInt16,
        gdal.GDT_UInt32,
        gdal.GDT_Float32,
        gdal.GDT_Float64,
    }

    # Read in data
    src_ds = gdal.Open(filename)
    srcband = src_ds.GetRasterBand(1)

    if no_data is None:
        no_data = srcband.GetNoDataValue()
        if no_data is None:
            raise Exception(
                "The source data did not have a NoData value. Please use the no_data argument to specify one. If should not be equal to any of the actual data values. If you are using all possible data values, then the situation is pretty hopeless - sorry."
            )

    srcdata = rdarray(srcband.ReadAsArray(), no_data=no_data)

    # raster_srs = osr.SpatialReference()
    # raster_srs.ImportFromWkt(raster.GetProjectionRef())

    if srcband.DataType not in allowed_types:
        raise Exception(
            "This datatype is not supported. Please file a bug report on RichDEM."
        )

    srcdata.projection = src_ds.GetProjectionRef()
    srcdata.geotransform = src_ds.GetGeoTransform()

    srcdata.metadata = dict()
    for k, v in src_ds.GetMetadata().items():
        srcdata.metadata[k] = v

    _AddAnalysis(
        srcdata, f"LoadGDAL(filename={filename}, no_data={no_data})"
    )

    return srcdata


def SaveGDAL(filename: str, rda: rdarray) -> None:
    """Save a GDAL file.

    Saves a RichDEM array to a data file in GeoTIFF format.

    If you need to do something more complicated, look at the source of this
    function.

    Args:
        filename (str):     Name of the raster file to be created
        rda      (rdarray): Data to save.

    Returns:
        No Return
    """
    if type(rda) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    if not GDAL_AVAILABLE:
        raise Exception("richdem.SaveGDAL() requires GDAL.")

    driver = gdal.GetDriverByName("GTiff")
    data_type = gdal.GDT_Float32  # TODO
    data_set = driver.Create(
        filename, xsize=rda.shape[1], ysize=rda.shape[0], bands=1, eType=data_type
    )
    data_set.SetGeoTransform(rda.geotransform)
    data_set.SetProjection(rda.projection)
    band = data_set.GetRasterBand(1)
    band.SetNoDataValue(rda.no_data)
    band.WriteArray(np.array(rda))
    for k, v in rda.metadata.items():
        data_set.SetMetadataItem(str(k), str(v))


def FillDepressions(dem: rdarray, epsilon: bool = False, in_place: bool = False, topology: str = "D8") -> Optional[rdarray]:
    """Fills all depressions in a DEM.

    Args:
        dem:       An elevation model
        epsilon:   If True, an epsilon gradient is imposed to all flat regions.
                     This ensures that there is always a local gradient.
        in_place:  If True, the DEM is modified in place and there is
                     no return; otherwise, a new, altered DEM is returned.
        topology:  A topology indicator

    Returns:
        DEM without depressions.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    if topology not in ["D8", "D4"]:
        raise Exception("Unknown topology!")

    if not in_place:
        dem = dem.copy()

    _AddAnalysis(dem, f"FillDepressions(dem, epsilon={epsilon})")

    demw = dem.wrap()

    if epsilon:
        if topology == "D8":
            _richdem.rdPFepsilonD8(demw)
        elif topology == "D4":
            _richdem.rdPFepsilonD4(demw)
    else:
        if topology == "D8":
            _richdem.rdFillDepressionsD8(demw)
        elif topology == "D4":
            _richdem.rdFillDepressionsD4(demw)

    dem.copyFromWrapped(demw)

    if not in_place:
        return dem


def BreachDepressions(dem: rdarray, in_place: bool = False, topology: str = "D8") -> Optional[rdarray]:
    """Breaches all depressions in a DEM.

    Args:
        dem     (rdarray): An elevation model
        in_place (bool):   If True, the DEM is modified in place and there is
                           no return; otherwise, a new, altered DEM is returned.
        topology (string): A topology indicator

    Returns:
        DEM without depressions.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    if topology not in ["D8", "D4"]:
        raise Exception("Unknown topology!")

    if not in_place:
        dem = dem.copy()

    _AddAnalysis(dem, "BreachDepressions(dem)")

    demw = dem.wrap()

    if topology == "D8":
        _richdem.rdBreachDepressionsD8(demw)
    elif topology == "D4":
        _richdem.rdBreachDepressionsD4(demw)

    dem.copyFromWrapped(demw)

    if not in_place:
        return dem


def ResolveFlats(dem: rdarray, in_place: bool = False) -> Optional[rdarray]:
    """Attempts to resolve flats by imposing a local gradient

    Args:
        dem          (rdarray):   An elevation model
        in_place (bool):   If True, the DEM is modified in place and there is
                           no return; otherwise, a new, altered DEM is returned.

    Returns:
        DEM modified such that all flats drain.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    if not in_place:
        dem = dem.copy()

    _AddAnalysis(dem, f"ResolveFlats(dem, in_place={in_place})")

    demw = dem.wrap()

    _richdem.rdResolveFlatsEpsilon(demw)

    dem.copyFromWrapped(demw)

    if not in_place:
        return dem


def FlowAccumulation(dem: rdarray, method: Optional[str] = None, exponent: Optional[float] = None, weights: Optional[rdarray] = None, in_place: bool = False) -> rdarray:
    """Calculates flow accumulation. A variety of methods are available.

    Args:
        dem      (rdarray): An elevation model
        method   (str):     Flow accumulation method to use. (See below.)
        exponent (float):   Some methods require an exponent; refer to the
                            relevant publications for details.
        weights  (rdarray): Flow accumulation weights to use. This is the
                            amount of flow generated by each cell. If this is
                            not provided, each cell will generate 1 unit of
                            flow.
        in_place (bool):    If True, then `weights` is modified in place. An
                            accumulation matrix is always returned, but it will
                            just be a view of the modified data if `in_place`
                            is True.

    =================== ============================== ===========================
    Method              Note                           Reference
    =================== ============================== ===========================
    Tarboton            Alias for Dinf.                `Taroboton (1997)              doi: 10.1029/96WR03137             <http://dx.doi.org/10.1029/96WR03137>`_
    Dinf                Alias for Tarboton.            `Taroboton (1997)              doi: 10.1029/96WR03137             <http://dx.doi.org/10.1029/96WR03137>`_
    Quinn               Holmgren with exponent=1.      `Quinn et al. (1991)           doi: 10.1002/hyp.3360050106        <http://dx.doi.org/10.1002/hyp.3360050106>`_
    Holmgren(E)         Generalization of Quinn.       `Holmgren (1994)               doi: 10.1002/hyp.3360080405        <http://dx.doi.org/10.1002/hyp.3360080405>`_
    Freeman(E)          TODO                           `Freeman (1991)                doi: 10.1016/0098-3004(91)90048-I  <http://dx.doi.org/10.1016/0098-3004(91)90048-I>`_
    FairfieldLeymarieD8 Alias for Rho8.                `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    FairfieldLeymarieD4 Alias for Rho4.                `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    Rho8                Alias for FairfieldLeymarieD8. `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    Rho4                Alias for FairfieldLeymarieD4. `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    OCallaghanD8        Alias for D8.                  `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    OCallaghanD4        Alias for D8.                  `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    D8                  Alias for OCallaghanD8.        `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    D4                  Alias for OCallaghanD4.        `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    =================== ============================== ===========================

    **Methods marked (E) require the exponent argument.**

    Returns:
        A flow accumulation according to the desired method. If `weights` was
        provided and `in_place` was True, then this matrix is a view of the
        modified data.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    facc_methods: Dict[str, Any] = {
        "Tarboton": _richdem.FA_Tarboton,
        "Dinf": _richdem.FA_Tarboton,
        "Quinn": _richdem.FA_Quinn,
        "FairfieldLeymarieD8": _richdem.FA_FairfieldLeymarieD8,
        "FairfieldLeymarieD4": _richdem.FA_FairfieldLeymarieD4,
        "Rho8": _richdem.FA_Rho8,
        "Rho4": _richdem.FA_Rho4,
        "OCallaghanD8": _richdem.FA_OCallaghanD8,
        "OCallaghanD4": _richdem.FA_OCallaghanD4,
        "D8": _richdem.FA_D8,
        "D4": _richdem.FA_D4,
    }

    facc_methods_exponent: Dict[str, Any] = {
        "Freeman": _richdem.FA_Freeman,
        "Holmgren": _richdem.FA_Holmgren,
    }

    if weights is not None and in_place:
        accum = rdarray(weights, no_data=-1)
    elif weights is not None and not in_place:
        accum = rdarray(weights, copy=True, meta_obj=dem, no_data=-1)
    elif weights is None:
        accum = rdarray(
            np.ones(shape=dem.shape, dtype="float64"), meta_obj=dem, no_data=-1
        )
    else:
        raise Exception("Execution should never reach this point!")

    if accum.dtype != "float64":
        raise Exception("Accumulation array must be of type 'float64'!")

    accumw = accum.wrap()

    _AddAnalysis(
        accum,
        "FlowAccumulation(dem, method={method}, exponent={exponent}, weights={weights}, in_place={in_place})".format(
            method=method,
            exponent=exponent,
            weights="None" if weights is None else "weights",
            in_place=in_place,
        ),
    )

    if method in facc_methods:
        facc_methods[method](dem.wrap(), accumw)
    elif method in facc_methods_exponent:
        if exponent is None:
            raise Exception(
                f'FlowAccumulation method "{method}" requires an exponent!'
            )
        facc_methods_exponent[method](dem.wrap(), accumw, exponent)
    else:
        raise Exception(
            "Invalid FlowAccumulation method. Valid methods are: "
            + ", ".join(list(facc_methods.keys()) + list(facc_methods_exponent.keys()))
        )

    accum.copyFromWrapped(accumw)

    return accum


def FlowAccumFromProps(props: rdarray, weights: Optional[rdarray] = None, in_place: bool = False) -> rdarray:
    """Calculates flow accumulation from flow proportions.

    Args:
        props    (rdarray): An elevation model
        weights  (rdarray): Flow accumulation weights to use. This is the
                            amount of flow generated by each cell. If this is
                            not provided, each cell will generate 1 unit of
                            flow.
        in_place (bool):    If True, then `weights` is modified in place. An
                            accumulation matrix is always returned, but it will
                            just be a view of the modified data if `in_place`
                            is True.

    Returns:
        A flow accumulation array. If `weights` was provided and `in_place` was
        True, then this matrix is a view of the modified data.
    """
    if type(props) is not rd3array:
        raise Exception("A richdem.rd3array or numpy.ndarray is required!")

    if weights is not None and in_place:
        accum = rdarray(weights, no_data=-1)
    elif weights is not None and not in_place:
        accum = rdarray(weights, copy=True, meta_obj=props, no_data=-1)
    elif weights is None:
        accum = rdarray(
            np.ones(shape=props.shape[0:2], dtype="float64"), meta_obj=props, no_data=-1
        )
    else:
        raise Exception("Execution should never reach this point!")

    if accum.dtype != "float64":
        raise Exception("Accumulation array must be of type 'float64'!")

    accumw = accum.wrap()

    _AddAnalysis(
        accum,
        "FlowAccumFromProps(dem, weights={weights}, in_place={in_place})".format(
            weights="None" if weights is None else "weights", in_place=in_place
        ),
    )

    _richdem.FlowAccumulation(props.wrap(), accumw)

    accum.copyFromWrapped(accumw)

    return accum


def FlowProportions(dem: rdarray, method: Optional[str] = None, exponent: Optional[float] = None) -> rdarray:
    """Calculates flow proportions. A variety of methods are available.

    Args:
        dem      (rdarray): An elevation model
        method   (str):     Flow accumulation method to use. (See below.)
        exponent (float):   Some methods require an exponent; refer to the
                            relevant publications for details.

    =================== ============================== ===========================
    Method              Note                           Reference
    =================== ============================== ===========================
    Tarboton            Alias for Dinf.                `Taroboton (1997)              doi: 10.1029/96WR03137             <http://dx.doi.org/10.1029/96WR03137>`_
    Dinf                Alias for Tarboton.            `Taroboton (1997)              doi: 10.1029/96WR03137             <http://dx.doi.org/10.1029/96WR03137>`_
    Quinn               Holmgren with exponent=1.      `Quinn et al. (1991)           doi: 10.1002/hyp.3360050106        <http://dx.doi.org/10.1002/hyp.3360050106>`_
    Holmgren(E)         Generalization of Quinn.       `Holmgren (1994)               doi: 10.1002/hyp.3360080405        <http://dx.doi.org/10.1002/hyp.3360080405>`_
    Freeman(E)          TODO                           `Freeman (1991)                doi: 10.1016/0098-3004(91)90048-I  <http://dx.doi.org/10.1016/0098-3004(91)90048-I>`_
    FairfieldLeymarieD8 Alias for Rho8.                `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    FairfieldLeymarieD4 Alias for Rho4.                `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    Rho8                Alias for FairfieldLeymarieD8. `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    Rho4                Alias for FairfieldLeymarieD4. `Fairfield and Leymarie (1991) doi: 10.1029/90WR02658             <http://dx.doi.org/10.1029/90WR02658>`_
    OCallaghanD8        Alias for D8.                  `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    OCallaghanD4        Alias for D8.                  `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    D8                  Alias for OCallaghanD8.        `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    D4                  Alias for OCallaghanD4.        `O'Callaghan and Mark (1984)   doi: 10.1016/S0734-189X(84)80011-0 <http://dx.doi.org/10.1016/S0734-189X(84)80011-0>`_
    =================== ============================== ===========================

    **Methods marked (E) require the exponent argument.**

    Returns:
        A flow proportion according to the desired method.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    fprop_methods = {
        "Tarboton": _richdem.FM_Tarboton,
        "Dinf": _richdem.FM_Tarboton,
        "Quinn": _richdem.FM_Quinn,
        "FairfieldLeymarieD8": _richdem.FM_FairfieldLeymarieD8,
        "FairfieldLeymarieD4": _richdem.FM_FairfieldLeymarieD4,
        "Rho8": _richdem.FM_Rho8,
        "Rho4": _richdem.FM_Rho4,
        "OCallaghanD8": _richdem.FM_OCallaghanD8,
        "OCallaghanD4": _richdem.FM_OCallaghanD4,
        "D8": _richdem.FM_D8,
        "D4": _richdem.FM_D4,
    }

    fprop_methods_exponent = {
        "Freeman": _richdem.FM_Freeman,
        "Holmgren": _richdem.FM_Holmgren,
    }

    fprops = rd3array(
        np.zeros(shape=dem.shape + (9,), dtype="float32"), meta_obj=dem, no_data=-2
    )
    fpropsw = fprops.wrap()

    _AddAnalysis(
        fprops,
        f"FlowProportions(dem, method={method}, exponent={exponent})",
    )

    if method in fprop_methods:
        fprop_methods[method](dem.wrap(), fpropsw)
    elif method in fprop_methods_exponent:
        if exponent is None:
            raise Exception(
                'FlowProportions method "' + method + '" requires an exponent!'
            )
        fprop_methods_exponent[method](dem.wrap(), fpropsw, exponent)
    else:
        raise Exception(
            "Invalid FlowProportions method. Valid methods are: "
            + ", ".join(
                list(fprop_methods.keys()) + list(fprop_methods_exponent.keys())
            )
        )

    fprops.copyFromWrapped(fpropsw)

    return fprops


def TerrainAttribute(dem: rdarray, attrib: str, zscale: float = 1.0) -> rdarray:
    """Calculates terrain attributes. A variety of methods are available.

    Args:
        dem    (rdarray):  An elevation model
        attrib (str):      Terrain attribute to calculate. (See below.)
        zscale (float):    How much to scale the z-axis by prior to calculation

    ======================= =========
    Method                  Reference
    ======================= =========
    slope_riserun           `Horn (1981)                   doi: 10.1109/PROC.1981.11918 <http://dx.doi.org/10.1109/PROC.1981.11918>`_
    slope_percentage        `Horn (1981)                   doi: 10.1109/PROC.1981.11918 <http://dx.doi.org/10.1109/PROC.1981.11918>`_
    slope_degrees           `Horn (1981)                   doi: 10.1109/PROC.1981.11918 <http://dx.doi.org/10.1109/PROC.1981.11918>`_
    slope_radians           `Horn (1981)                   doi: 10.1109/PROC.1981.11918 <http://dx.doi.org/10.1109/PROC.1981.11918>`_
    aspect                  `Horn (1981)                   doi: 10.1109/PROC.1981.11918 <http://dx.doi.org/10.1109/PROC.1981.11918>`_
    curvature               `Zevenbergen and Thorne (1987) doi: 10.1002/esp.3290120107  <http://dx.doi.org/10.1002/esp.3290120107>`_
    planform_curvature      `Zevenbergen and Thorne (1987) doi: 10.1002/esp.3290120107  <http://dx.doi.org/10.1002/esp.3290120107>`_
    profile_curvature       `Zevenbergen and Thorne (1987) doi: 10.1002/esp.3290120107  <http://dx.doi.org/10.1002/esp.3290120107>`_
    ======================= =========

    Returns:
        A raster of the indicated terrain attributes.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    terrain_attribs = {
        # "spi":                _richdem.TA_SPI,
        # "cti":                _richdem.TA_CTI,
        "slope_riserun": _richdem.TA_slope_riserun,
        "slope_percentage": _richdem.TA_slope_percentage,
        "slope_degrees": _richdem.TA_slope_degrees,
        "slope_radians": _richdem.TA_slope_radians,
        "aspect": _richdem.TA_aspect,
        "curvature": _richdem.TA_curvature,
        "planform_curvature": _richdem.TA_planform_curvature,
        "profile_curvature": _richdem.TA_profile_curvature,
    }

    if attrib not in terrain_attribs:
        raise Exception(
            "Invalid TerrainAttributes attribute. Valid attributes are: "
            + ", ".join(terrain_attribs.keys())
        )

    result = rdarray(
        np.zeros(shape=dem.shape, dtype="float32"), meta_obj=dem, no_data=-9999
    )
    resultw = result.wrap()

    _AddAnalysis(
        result,
        f"TerrainAttribute(dem, attrib={attrib}, zscale={zscale})"
    )

    terrain_attribs[attrib](dem.wrap(), resultw, zscale)

    result.copyFromWrapped(resultw)

    return result

def generate_perlin_terrain(size: int, seed: int) -> rdarray:
    """Generates random terrain based on Perlin noise

    Args:
        size: Size of one edge of the terrain. The terrain will be a square.
        seed: Random seed used to generate the terrain.

    Returns:
        A DEM
    """
    dem = rdarray(np.zeros(shape=(size, size), dtype="float64"), no_data=-9999)
    dem.geotransform = STANDARD_GEOTRANSFORM
    demw = dem.wrap()
    _richdem.generate_perlin_terrain(demw, seed)
    _AddAnalysis(dem, f"GeneratePerlinTerrain(size={size}, seed={seed})")
    dem.copyFromWrapped(demw)
    return dem

def get_depression_hierarchy(dem: rdarray, labels: rdarray) -> Tuple[List[depression_hierarchy.Depression], rdarray]:
    """Fills all depressions in a DEM.

    Args:
        dem:       An elevation model
        labels:    Should have `OCEAN` for cells representing the "ocean" (the
                   place to which depressions drain) and `NO_DEP` for all other
                   cells.

    Returns:
        label:  Modified in-place to contain a label indicating which depression
                the cell belongs to. The indicated label is always the leaf of
                the depression hierarchy, or the OCEAN.

        flowdirs:  A value [0,7] indicated which direction water from the cell
                   flows in order to go "downhill". All cells have a flow
                   direction (even flats) except for pit cells.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required!")

    demw = dem.wrap()
    labelsw = labels.wrap()

    flowdirs = rdarray(_richdem.NO_FLOW * np.ones(dem.shape, dtype=np.int8), no_data=-9999, geotransform=STANDARD_GEOTRANSFORM)
    flowdirsw = flowdirs.wrap()

    dhret = depression_hierarchy.get_depression_hierarchy(demw, labelsw, flowdirsw)

    return dhret, flowdirs

def get_new_depression_hierarchy_labels(shape: Tuple[int, int], no_data: float = -9999, geotransform: Optional[np.ndarray] = None) -> rdarray:
    """Get a new labels array with for use with the depression hierarchy

    Will have ocean on the borders and NO_DEP everywhere else.

    Args:
        no_data: `no_data` value
        geotransform: A GDAL geotransform; defaults to `richdem.STANDARD_GEOTRANSFORM`.

    Returns:
        The new labels array
    """
    dhlret = rdarray(
      depression_hierarchy.OCEAN * np.ones(shape, dtype = np.uint32),
      no_data=no_data,
      geotransform=geotransform if geotransform is not None else STANDARD_GEOTRANSFORM
    )
    dhlret[1:-1,1:-1] = depression_hierarchy.NO_DEP
    return dhlret

def fill_spill_merge(dem: rdarray, labels: rdarray, flowdirs: rdarray, deps: List[depression_hierarchy.Depression], wtd: rdarray) -> None:
    """
    This function routes surface water into pit cells and then distributes
    it so that it fills the bottoms of depressions, taking account of overflows.

    Args:
    topo     Topography used to generate the DepressionHierarchy
    label    Labels from GetDepressionHierarchy indicate which depression
             each cell belongs to.
    flowdirs Flowdirs generated by GetDepressionHierarchy
    deps     The DepressionHierarchy generated by GetDepressionHierarchy
    wtd      Water table depth. Values of 0 indicate saturation.
             Negative values indicate additional water can be added to the
             cell. Positive values indicate standing surface water.

    Note that from get_depression_hierarchy we already know the number of cells
    within each depression as well as its volume.

    @return Modifies the depression hierarchy `deps` to indicate the amount of
            water contained in each depression. Modifies `wtd` to indicate how
            saturated a cell is or how much standing surface water it has.
    """
    if type(dem) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required for the argument `dem`!")
    if type(labels) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required for the argument `labels`!")
    if type(flowdirs) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required for the argument `flowdirs`!")
    if type(wtd) is not rdarray:
        raise Exception("A richdem.rdarray or numpy.ndarray is required for the argument `wtd`!")

    demw = dem.wrap()
    labelsw = labels.wrap()
    flowdirsw = flowdirs.wrap()
    wtdw = wtd.wrap()

    dhret = depression_hierarchy.fill_spill_merge(demw, labelsw, flowdirsw, deps, wtdw)
