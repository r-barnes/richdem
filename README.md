RichDem
=======

Author: Richard Barnes (rbarnes@umn.edu)

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly
process even very large DEMs.

**Note: Further development is occurring on another branch. RichDEM will shortly be rereleased as a suite of easy-to-use libraries and utility programs.**

RichDEM can use both the D8 and D-infinite (Tarboton) flow metrics. It can
resolve terrain depressions (or pits) either by filling or by channel carving.
It can calculate contributing/up-slope areas, slopes, curvatures, and aspects.

Future versions may use Intel's Thread Building Blocks to achieve additional
increases in speed.

Using It
========

Casual users will be interested in the **dist/** directory which contains
pre-compiled Windows binaries and an ArcGIS toolbox for performing pit-filling,
flow accumulation, CTI, and SPI caclulations under both the D8 and D-finite
paradigms.

More advanced users will appreciate the ability to run **DOxygen** on
**doxy.conf** to generate thorough documentation. Browsing "main.cpp" should be
sufficient to familiarize oneself with many of the functions and interface
styles of the code. Throughout, the code is written in what I hope is an
easy-to-follow style and many functions include extensive (and growing!)
documentation.

Specific Algorithms
===================
Many of the algorithms used in RichDEM are documented in journal or conference
publications. In the case of older algorithms by other authors, it is often
possible to find the paper in the literature. Some of the newer algorithms I
developed have not yet had a chance to be widely utilized. These algorithms are
listed below.

Additionally, each publication has its own GitHub repository featuring
easily-compiled, minimum working examples of the algorithms, along with examples
and test data sets.

These are available as follows:

 * Flat-resolution algorithm. [Link](https://github.com/r-barnes/Barnes2013-FlatSurfaces)
 * Depression-filling algorithm. [Link](https://github.com/r-barnes/Barnes2013-Depressions)

Publications
============
The algorithms used in RichDEM have been published in the following articles:

* Barnes, Lehman, Mulla. "An Efficient Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.01.009](http://dx.doi.org/10.1016/j.cageo.2013.01.009)

* Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.04.024](http://dx.doi.org/10.1016/j.cageo.2013.04.024)

* Barnes, Lehman, Mulla. "Distributed Parallel D8 Up-Slope Area Calculation in Digital Elevation Models". Intn'l Conf. on Parallel & Distributed Processing Techniques & Applications. [Link](http://rbarnes.org/section/sci/2011_barnes_distributed.pdf)

Requires:
 * OpenMP
