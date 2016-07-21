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



Design Philosophy
=================

The design of RichDEM is guided by these principles:

* **Algorithms will be well-tested.** Every algorithm is verified by a rigorous
  testing procedure. See below.

* **Algorithms will be fast, without compromising safety and accuracy.** The
  algorithms used in RichDEM are state of the art, permitting analyses that
  would take days on other systems to be performed in hours, or even minutes.

* **Algorithms will be available as libraries, whenever possible.** RichDEM is
  designed as a set of header-only C++ libraries, making it easy to include in
  your projects and easy to incorporate into other programming languages.
  RichDEM also includes apps, which are simple wrappers around the algorithms, 
  and a limited, but growing, set of algorithms which may have special
  requirements, like MPI, that make them unsuitable as libraries. These are 
  available as programs.

* **Programs will have a command-line interface, not a GUI.** Command-line
  interfaces are simple to use and offer extreme flexibility for both users and
  programmers. They are available on every type of operating system. RichDEM
  does not officially support any GUI. Per the above, encapsulating RichDEM in
  a high-level interface of your own is not difficult.

* **Algorithms and programs will be portable.** Linux, Mac, and Windows should
  all be supported.

* **The code will be beautiful.** RichDEM's code utilizes sensible variable
  names and reasonable abstractions to make it easy to understand, use, and
  design algorithms. The code contains extensive internal documentation which is
  DOxygen compatible.

* **Programs and algorithms will provide useful feedback.** Progress bars will 
  appear if desired and the output will be optimized for machine parsing.



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

Testing Methodology
===================

TODO

Parsable Output
===================

Every line of output from RichDEM begins with one of the following characters,
making it easy to parse with a machine.

 * **c**: Configuration information: program version, input files, and command
          line options, &c.

 * **p**: Progress information: inform the user to keep calm because we're
          carrying on.

 * **t**: Timing information: How long stuff took

 * **i**: I/O: Amount of data loaded from disk

 * **n**: I/O: Amount of data transferred through a network

 * **r**: Amount of RAM used

 * **A**: Algorithm name

 * **C**: Citation for algorithm

 * **m**: Miscallaneous counts

 * **E**: Indicates an error condition

 * **W**: Indicates a warning

 * **d**: Debugging info

All output data shall have the form:

    <INDICATOR CHARACTER> <MESSAGE/MEASUREMENT NAME> [= <VALUE> [UNIT]]

The amount of whitespce may very for aesthetic purposes.

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

* Barnes, R., 2016. Parallel priority-flood depression filling for trillion cell digital elevation models on desktops or clusters. Computers & Geosciences. doi:[10.1016/j.cageo.2016.07.001](http://dx.doi.org/10.1016/j.cageo.2016.07.001)

* Barnes, Lehman, Mulla. "An Efficient Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.01.009](http://dx.doi.org/10.1016/j.cageo.2013.01.009)

* Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models". Computers &amp; Geosciences, 2013. doi: [10.1016/j.cageo.2013.04.024](http://dx.doi.org/10.1016/j.cageo.2013.04.024)

* Barnes, Lehman, Mulla. "Distributed Parallel D8 Up-Slope Area Calculation in Digital Elevation Models". Intn'l Conf. on Parallel & Distributed Processing Techniques & Applications. [Link](http://rbarnes.org/section/sci/2011_barnes_distributed.pdf)

Requires:
 * OpenMP
