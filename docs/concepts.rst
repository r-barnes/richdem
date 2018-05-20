Concepts
===================================

Gridded Data
-----------------------------------

RichDEM assumes that data is provided in the form of a rectangular grid of cells
with some *width* and *height*. Furthermore, the data comprising this grid must
be laid out in a flat array such that the value of any cell `(x,y)` can be
accessed via the equation `y*width+x`. This is known as **row-major** ordering.

Data can be passed to RichDEM in one of three ways.

  1. Data can be loaded via GDAL. GDAL handles the heavy lifting of ensuring
     that data is loaded in a form which complies with the above assumptions.

  2. Data can be manually added to a `richdem::Array2D<T>` object (C++) or a 
     `richdem.rdarray` object (Python).

  3. A `richdem::Array2D<T>` (C++) or `richdem.rdarray` (Python) object can be
     used to wrap existing row-major memory. This capability allows RichDEM to
     easily integrate with existing code.



Metadata
-----------------------------------

A RichDEM array is accompanied by several kinds of metadata. These can be loaded
by GDAL or specified manually.

 - A **NoData value**, as discussed below.
 - A **projection**. This is, typically, a PROJ4 or WKT string that identifies 
   the projection the data maps to. The choice of projection does not affect
   RichDEM's operations.
 - A **geotransform**. This is a six element array which determines where in a
   projection a RichDEM array's data is located, as well as the cell sizes.
   Further details are below. This setting does affect how RichDEM processes
   data.
 - A **metadata** entry which contains arbitrary metadata strings such as 
   `PROCESSING_HISTORY` (see below).



Geotransform
-----------------------------------

A **geotransform** is a six element array which determines where in a
projection a RichDEM array's data is located, as well as the cell sizes.

Typically, the geotransform is an affine transform consisting of six
coefficients which map pixel/line coordinates into a georeferenced space using
the following relationship:

    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)

In case of north up images, the `GT(2)` and `GT(4)` coefficients are zero, and
the `GT(1)` is pixel width, and `GT(5)` is pixel height. The `(GT(0),GT(3))`
position is the top left corner of the top left pixel of the raster.

Note that the pixel/line coordinates in the above are from `(0.0,0.0)` at the
top left corner of the top left pixel to (width_in_pixels,height_in_pixels) at
the bottom right corner of the bottom right pixel. The pixel/line location of
the center of the top left pixel would therefore be `(0.5,0.5)`

(Text drawn from GDAL documentation.)



NoData values
-----------------------------------

RichDEM recognizes cells with a *NoData* value and treats them in special ways.
The *NoData* value is a number, such as `-9999` which represents a cell that
isn't part of a data set.

In depression-filling and the determination of flow directions, *NoData* cells
are treated as having a lower elevation than any other cell, which enforces
drainage to the edge of the DEM.

Because *NoData* values fundamental affect operations, RichDEM requires that you
specify what *NoData* value it should use.

It is important that you choose a value that doesn't correspond to any of your
actual data values. In 1-byte/8-bit DEMs this may not be possible since there
are only 256 distinct values available. In this case, you should cast your data
to a 2-byte/16-bit form and choose a *NoData* value that is not in the range
`0-255`.

Other terrain analysis programs may use a binary masking array to indicate which
cells are *NoData*. RichDEM does not do this because usually a minority of cells
are *NoData* and a separate array increases the amount of memory used and
reduces cache utilization, both of which may reduce performance.

================= =================== ========================
Language          Set                 Get
================= =================== ========================
Python            `rda.no_data=-1`    `rda.no_data`
C++               `rda.setNoData(-1)` `rda.noData( )`
================= =================== ========================



Processing History
-----------------------------------

RichDEM automagically keeps track of what operations are performed on your data.
This means that your outputs will contain an exact record of what operations
were used to produce them. This helps ensure reproducibility. And, remember,
good science is reproducible science.

In the following, a series of operations is performed and the Processing History
is then examined.

.. code-block:: python

    import richdem as rd
    dem   = rd.LoadGDAL("../data/beauford.tif")
    rd.FillDepressions(dem, epsilon=False, in_place=True)
    accum = rd.FlowAccumulation(dem, method='D8')
    rd.SaveGDAL('/z/out.tif', accum)
    print(accum.metadata['PROCESSING_HISTORY'])

The processing history of a saved dataset can be viewed using a few different
commands:

.. code-block:: bash

    gdalinfo /z/out.tif
    rd_info  /z/out.tif

The processing history appears as follows.

.. code-block:: text

    2017-12-20 17:55:19.892388 UTC | RichDEM (Python 0.0.4) (hash=e02d5e2, hashdate=2017-12-19 23:52:52 -0600) | LoadGDAL(filename=../data/beauford.tif, no_data=-9999.0)
    2017-12-20 17:55:19.900234 UTC | RichDEM (Python 0.0.4) (hash=e02d5e2, hashdate=2017-12-19 23:52:52 -0600) | FillDepressions(dem, epsilon=False)
    2017-12-20 17:55:20.514098 UTC | RichDEM (Python 0.0.4) (hash=e02d5e2, hashdate=2017-12-19 23:52:52 -0600) | FlowAccumulation(dem, method=D8)

Note that the **first column** is the time at which the operation was performed,
the **second column** is the program which performed the operation, and the
**third column** is the command which was run.

================= ==============================
Language          Command
================= ==============================
Python            `rda.metadata`
C++               `rda.metadata`
================= ==============================


In-Place Operations
-----------------------------------

To save memory RichDEM performs some operations, such as depression-filling, in
place. This means that the data is modified and the original data will be lost
unless it has been copied.

For instance, in Python `FillDepressions` has two distinct forms:

.. code-block:: python

    #In-place filling, no return value
    rd.FillDepressions(dem, in_place=True)
    #Fill a copy
    dem_filled = rd.FillDepressions(dem, in_place=False)

whereas in C++, a copy must be made:

.. code-block:: python

    #In-place filling, no return value
    richdem::FillDepressions(dem)

    #Fill a copy
    auto demcopy = dem; //TODO: Make sure this syntax is right
    richdem::FillDepressions(demcopy)


Topology
-----------------------------------

RichDEM offers two topologies, though not all functions differentiate between
them. Thse are:

 - **D8**: The cells are arranged in a regular, rectilinear grid. Each cell
           connects with each of its neighbouring cells.

 - **D4**: The cells are arranged in a regular, rectilnear grid. Each cell
           connects with the cells to its north, south, east, and west (the
           cells up, down, left, and right of it).

In C++, the foregoing topologies are accessed via the :code:`Topology`
enumeration, similar to the following:

.. code-block:: c++

    FillDepressions<Topology::D8>(dem);
    FillDepressions<Topology::D4>(dem);
