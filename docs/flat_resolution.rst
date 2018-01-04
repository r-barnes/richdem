Flat Resolution
===============

The problem with doing complete filling on flats is that the resulting DEM
contains mathematically flat areas with no local gradient. This makes it
impossible to determine flow directions for these areas.

Imposing an epsilon gradient during depression-filling is one solution as
discussed in Epsilon Filling (:ref:`epsilon-filling-label`); however, the
gradient produced may appear unaesthetic because drainage takes a least-distance
route to the flat's edges.

We stress the word **aesthetic** here since, after depression-filling, no
information about local gradients remains from the original DEM, so, in a sense,
all reconstructed drainage patterns are equally silly. Sometimes, though, this
is the best you can do.

This pages discusses alternatives.



Barnes (2014) Flat Resolution
-----------------------------

    Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:10.1016/j.cageo.2013.01.009

The Barnes (2014) flat resolution algorithm routes flow both towards the edges
of flats and away from the high areas surrounding them. The result is an
aesthetically pleasing drainge pattern.

It can do so either by adjust the elevations of the DEM's cells or by adjusting
flow metrics derived from the DEM.



Elevation Adjustment
~~~~~~~~~~~~~~~~~~~~

In this method, the elevation of a DEM can be adjusted so that every cell in the
region is raised some small amount, ε, above cells which are closer to a
depression's spill point and farther from its surrounding high areas.

This must be done carefully. In floating-point DEMs, the value ε is non-constant
and must be chosen using the !`std::nextafter` function. If a depression is too
large, the imposed gradient may result in the interior of the depression being
raised above the surrounding landscape. Using `double` instead of `float`
reduces the potential for problems at a cost of twice the space used. If a
problem does arise, RichDEM provides a warning.

Recall from Epsilon Filling (:ref:`epsilon-filling-label`) that an epsilon
gradient imposed during depression-filling results in an elevation adjustment
that looks like this:

.. plot::
    :width: 800pt
    :include-source:
    :context: reset
    :outname: flat_resolution_dep_fill_epsilon

    #Load dataset
    beau        = rd.rdarray(np.load('imgs/beauford.npz')['beauford'], no_data=-9999)

    #Fill the depression entirely
    beau_filled = rd.FillDepressions(beau, epsilon=False, in_place=False)

    #Construct the epsilon drainage surface via filling
    beau_eps    = rd.FillDepressions(beau, epsilon=True, in_place=False)

    diff    = beau_eps - beau
    rd.rdShow(diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(8,5.5))

In contrast, the Barnes (2014) convergent elevation adjustment looks like this:

.. plot::
    :width: 800pt
    :include-source:
    :context: close-figs
    :outname: flat_resolution_barnes2014_epsilon

    #Resolve flats by imposing a convergent epsilon gradient
    beau_flat_eps = rd.ResolveFlats(beau_filled, in_place=False)

    diff    = beau_flat_eps - beau
    rd.rdShow(diff, ignore_colours=[0], axes=False, cmap='jet', figsize=(8,5.5))

================= ============================================
Language          Command
================= ============================================
Python            `richdem.rdResolveFlats()`
C++               `richdem::ResolveFlatsEpsilon()`
================= ============================================

+-------------------+--------------------------------------------+
|Pros               | Cons                                       |
+-------------------+--------------------------------------------+
| - All cells drain | - Not as fast as simple depression filling |
|                   | - May modify large portions of a DEM       |
|                   | - May create elevated regions              |
|                   | - Success may depend on data type          |
+-------------------+--------------------------------------------+



Flow Metric Adjustment
----------------------

.. todo:: TODO: Flow metric adjustment

TODO