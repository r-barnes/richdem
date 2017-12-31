Example DEMs
================================

Beauford Watershed, Minnesota, USA is frqeuently used as an example dataset.

.. plot::
    :width: 800pt
    :include-source:

    import richdem as rd
    beau    = rd.LoadGDAL("/home/rick/data/gis/beauford.tif")
    rd.rdShow(beau, ignore_colours=[0], axes=False, cmap='jet', figsize=(8,5.5))
