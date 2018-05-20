Depression-Breaching
====================

Depressions, otherwise known as pits, are areas of a landscape wherein flow
ultimately terminates without reaching an ocean or the edge of a digital
elevation model.



Depressions, Pits, and Sinks
----------------------------

Depressions have been called by a variety of names. To clarify this mess,
Lindsay (2016) provides a typology. This typology is followed here.

.. image:: imgs/Lindsay2016_depression_typology.png


Original DEM
----------------------------

For reference, the original DEM appears as follows:

.. plot::
    :width: 800pt
    :include-source:
    :context: reset
    :outname: depression_breach_original

    import richdem as rd
    import numpy as np
    
    beau    = rd.rdarray(np.load('imgs/beauford.npz')['beauford'], no_data=-9999)
    beaufig = rd.rdShow(beau, ignore_colours=[0], axes=False, cmap='jet', figsize=(8,5.5))


Complete Breaching
----------------------------

Depression-breaching is used to dig channels from the pit cells of a DEM to the
nearest cells (in priority-flood sense) outside of the depression containing the
pit. This resolves the depression as all cells in the depression now have a
drainage path to the edge of the DEM.

The result looks as follows:

.. plot::
    :width: 800pt
    :include-source:
    :context: close-figs
    :outname: breaching_complete

    beau_breached    = rd.BreachDepressions(beau, in_place=False)
    beaufig_breached = rd.rdShow(beau_breached, ignore_colours=[0], axes=False, cmap='jet', vmin=beaufig['vmin'], vmax=beaufig['vmax'], figsize=(8,5.5))

We can visualize the difference between the two like so:

.. plot::
    :width: 800pt
    :include-source:
    :context: close-figs
    :outname: depression_complete_breached_original_diff

    beau_diff    = beau_breached - beau
    beaufig_diff = rd.rdShow(beau_diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(8,5.5))

Complete Breaching is available via the following commands:

================= ==============================
Language          Command
================= ==============================
Python            `richdem.BreachDepressions()`
C++               `richdem::BreachDepressions<Topology>()`
================= ==============================

+-------------------------------+--------------------------------------+
|Pros                           |  Cons                                |
+-------------------------------+--------------------------------------+
| - Minimal Modifcations to DEM | - Slightly slower                    |
| - Simple                      |                                      |
+-------------------------------+--------------------------------------+

