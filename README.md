RichDem
=======

Author: Richard Barnes (rbarnes@umn.edu)

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly
process even very large DEMs.

RichDEM can use both the D8 and D-infinite (Tarboton) flow metrics. It can
resolve terrain depressions (or pits) either by filling or by channel carving.
It can calculate contributing/up-slope areas, slopes, curvatures, and aspects.

Future versions may use Intel's Thread Building Blocks to achieve additional
increases in speed.

To use RichDEM, edit "main.cpp" to suit your needs and then simply run "make".

To generate DOxygen documentation use "doxygen doxy.conf".

Publications
============
The algorithms used in RichDEM have been published in the following articles:

* Barnes, Lehman, Mulla. "An Efficient Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.01.009](http://dx.doi.org/10.1016/j.cageo.2013.01.009)

* Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.04.024](http://dx.doi.org/10.1016/j.cageo.2013.04.024)

* Barnes, Lehman, Mulla. "Distributed Parallel D8 Up-Slope Area Calculation in Digital Elevation Models". Intn'l Conf. on Parallel & Distributed Processing Techniques & Applications. [Link](http://rbarnes.org/section/sci/2011_barnes_distributed.pdf)

Requires:
 * OpenMP
