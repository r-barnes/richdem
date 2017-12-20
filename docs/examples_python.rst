Python Examples
=============================

Depression-filling a DEM and saving it
--------------------------------------

.. code-block:: python

    import richdem as rd

    #Load DEM
    dem     = rd.LoadGDAL("mydem.tif")
    
    #Fill depressions in the DEM. The data is modified in-place to avoid making
    #an unnecessary copy. This saves both time and RAM. Note that the function
    #has no return value when `in_place=True`.
    rd.FillDepressions(dem, epsilon=False, in_place=True)
    
    #Save the DEM
    rd.SaveGDAL("mydem_filled.tif", dem)



Comparing filled vs. unfilled DEMs
--------------------------------------

.. code-block:: python

    import richdem as rd

    #Load DEM
    dem     = rd.LoadGDAL("mydem.tif")
    
    #Copy the DEM so we can compare the altered DEM to the unaltered original
    demorig = dem.copy()
    
    #Fill depressions in the DEM. The data is modified in place, but, since we
    #made a copy above neither time nor memory is really saved.
    rd.FillDepressions(dem, epsilon=False, in_place=True)
    
    #Get the difference of the filled and unfilled DEM
    diff = dem - demorig
    
    #Display the difference. Do not plot values where there was no difference.
    rd.rdShow(diff, ignore_colours=[0])

The foregoing could be written more succinctly by using `in_place=False`:

.. code-block:: python

    import richdem as rd

    #Load DEM
    demorig = rd.LoadGDAL("mydem.tif")
    
    #Fill depressions in the DEM. By using `in_place=False`, a copy of the DEM
    #is made and only this copy is altered. Note that the function now has a 
    #return value.
    dem = rd.FillDepressions(dem, epsilon=False, in_place=False)
    


The rdarray class
--------------------------------------

RichDEM has a class `rdarray` which can wrap around NumPy arrays without copying
the memory therein. This makes it easy to pass data to RichDEM in a format it
understands. 

The `rdarray` class works exactly like a NumPy array, but has some special
features.

Namely, an `rdarray` has the following properties:

============== ================================================================
metadata       A dictionary containing metadata in key-value pairs.
projection     A string describing the geographic projection of the dataset
geotransform   An array of six floating-point values representing the size and offset of the DEM's cells.
no_data        A value indicating which cell values should be treated as special NoData cells. (See Concepts, TODO)
============== ================================================================

The `metadata` dictionary contains the special entry
`metadata['PROCESSING_HISTORY']`. This entry contains a complete list of
everything that RichDEM has done to a dataset. (See Concepts, TODO)



Using RichDEM without GDAL
--------------------------------------

GDAL is an optional dependency to RichDEM. In modeling, data is often stored in
NumPy arrays and evolved or manipulated without ever being loaded from or saved
to disk. To use NumPy arrays, simply wrap them with an `rdarray` as shown below.
Details about the `rdarray` are above.

Since an `rdarray` must have a `no_data` value and choosing the `no_data` value
incorrectly can produce erroneous results, RichDEM does not automagically
convert NumPy arrays. It must be done manually.

.. code-block:: python

    import richdem as rd
    import numpy as np

    #Create some NumPy data
    npa = np.random.random(size=(100,100))

    #Wrap the NumPy data in an rdarray. I want to treat all of the cells as data
    #cells, so I use `no_data=-9999` since I know none of my cells will have
    #this value.
    rda = rd.rdarray(npa, no_data=-9999)

    #Fill depressions, modifying in place. At this point, the calculation I
    #wanted to do is done and I can throw away the `rda` object.
    rd.FillDepressions(rda, in_place=True)