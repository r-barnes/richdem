Loading Data
========================================================

Python
--------------------------------------------------------

Data can be loaded in several ways.

To load from disk, if GDAL is available on your system, almost any form of
raster data can be easily loaded, like so:


GDAL
~~~~

.. code-block:: python

    import richdem as rd
    beau = rd.LoadGDAL("beauford.tif")


NumPy
~~~~~

Data can also be loaded from a NumPy array:

.. code-block:: python

    import numpy as np
    import richdem as rd

    npa = np.random.random(size=(50,50))
    rda = rd.rdarray(npa, no_data=-9999)

**Note** that !`rd.rdarray()` creates a *view* of the data stored in !`npa`.
Modifying `rda` will modify `npa`. This prevents unwanted memory from being
used. If you instead want `rda` to be a new copy of the data, use:

.. code-block:: python

    rda = rd.rdarray(a, no_data=-9999)


Saved NumPy Arrays
~~~~~~~~~~~~~~~~~~

It is possible to save, and load, data to and from a NumPy array like so:

.. code-block:: python

    import numpy as np
    import richdem as rd

    npa = np.random.random(size=(50,50))
    rda = rd.rdarray(npa, no_data=-9999)
    np.save('out.npy', rda)
    loaded = rd.rdarray(np.load('out.npy'), no_data=-9999)

This can be done in a compressed format like so:

.. code-block:: python

    np.savez('rda', rda=rda)
    np.load('rda.npz')['rda']

Note that there is not yet a way to save the metadata of an rdarray. (TODO)



C++
--------------------------------------------------------

TODO