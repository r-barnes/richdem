RichDEM
=======

[![ReadTheDocs](https://readthedocs.org/projects/richdem/badge/?version=latest)](https://richdem.readthedocs.io/)
[![Travis](https://travis-ci.org/r-barnes/richdem.svg?branch=master)](https://travis-ci.org/r-barnes/richdem)
[![DOI](https://zenodo.org/badge/7469158.svg)](https://zenodo.org/badge/latestdoi/7469158)

Author: Richard Barnes (rbarnes@umn.edu)

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly
process even very large DEMs.

RichDEM offers a variety of flow metrics, such as D8 and D∞. It can flood or
breach depressions. It can calculate flow accumulation, slops, curvatures, &c.

RichDEM is available as a performant C++ library, a low-dependency Python
package, and a set of command-line tools.

Please cite RichDEM (see below).



Using It
========

Citing It
---------

As of 883ea734e957, David A. Wheeler's SLOCCount estimates the value of RichDEM
at $240,481 and 1.78 person-years of development effort. This value is yours to
use, but citations are encouraged as they provide justification of continued
development.

General usage of the library can be cited as:

    Barnes, Richard. 2016. RichDEM: Terrain Analysis Software. http://github.com/r-barnes/richdem

An example BibTeX entry is:

    @manual{RichDEM,
      title        = {RichDEM: Terrain Analysis Software},
      author       = {Richard Barnes},
      year         = {2016},
      url          = {http://github.com/r-barnes/richdem}, 
    }

This information will be updated as versioned releases become available.

Although I have written all of the code in this library, some of the algorithms
were discovered or invented by others, and they deserve credit for their good
work. Citations to particular algorithms will be printed whenever an app,
program, or library function is run. Such citations are prefixed by the
character `C` and look like:

    C Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024

A typical academic citation might read as follows:

 > We performed hydrological corrections on our DEM using the Zhou (2016) algorithm implemented in RichDEM (Barnes 2016).

As A Python Package
-------------------

Get the package with:

    pip3 install richdem

And use:

    import richdem

The command:

    help(richdem)

provides all the relevant documentation.



As A Command-line Tool
----------------------

To get the command-line tools, install the Python package with:

    pip3 install richdem

The command-line tools are all named `rd_*`, so typing `rd_` on your command-
line and hitting tab a few times should give you the full list of what's
available.



As A Library
------------

Upon compilation, point your library search path to the `include` directory.
Include various files using, e.g.

    #include "richdem/common/Array2D.hpp"

All files include extensive documentation. At this stage the location of certain
functions may be subject to change. This will be noted in the `NEWS` file. (TODO)



As A Handy Collection of Tools
------------------------------

Running `make` in the `apps` directory will produce a large number of useful
scripts which are essentially wrappers around standard uses of the RichDEM
libraries. The [apps/README.md](apps/README.md) file and the apps themselves
contain documentation explaining what they all do.



For Processing Large Datasets
-----------------------------

The `programs` directory contains several programs which have not been converted
to libraries. This is usually because their functionality is specific and they
are unlikely to be useful as a library. Each directory contains a makefile and a
readme explaining the purpose of the program.



Documentation
-------------

Documentation is available at [richdem.com](http://richdem.com/doc/index.html).
The documentation is auto-generated from the many `README.md` files throughout
the codebase and the extensive comments in the source code.



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

* **Analyses will be reproducible.** Every time you run a RichDEM command that
  command is logged and timestamped in the output data, along with the version
  of the program you created the output with. Additionally, a history of all
  previous manipulations to the data is kept. Use `rd_view_processing_history`
  to see this.**



Testing Methodology
===================

Simple algorithms are shown to be correct through visual inspection and
comparison against hand-crafted examples. Correctness for more complex
algorithms is often "boot-strapped" by comparing the results of simple
algorithms versus the complex algorithms on a large number of randomly-generated
datasets.

This is a work in progress. TODO


Correctness
===========

Correctness is established via a number of methodologies building from code
inspection in the simplest cases to output comparison between simple and complex
implementations.

Correctness is noted in source code comments under `@correctness` sections.
These are, in turn, printed to the Doxygen documentation output.

A master list of how correctness was established for each algorithm is available
at [tests/README.md](tests/README.md).



Parsable Output
===================

Every line of output from RichDEM begins with one of the following characters,
making it easy to parse with a machine.

 * **A**: Algorithm name

 * **a**: Analysis command: the command line used to run the program

 * **c**: Configuration information: program version, input files, and command
          line options, &c.

 * **C**: Citation for algorithm

 * **d**: Debugging info

 * **E**: Indicates an error condition
 
 * **i**: I/O: Amount of data loaded from disk

 * **m**: Miscallaneous counts

 * **n**: I/O: Amount of data transferred through a network

 * **p**: Progress information: inform the user to keep calm because we're
          carrying on.
 
 * **r**: Amount of RAM used

 * **t**: Timing information: How long stuff took

 * **W**: Indicates a warning


All output data shall have the form:

    <INDICATOR CHARACTER> <MESSAGE/MEASUREMENT NAME> [= <VALUE> [UNIT]]

The amount of whitespace may very for aesthetic purposes.



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
 * Large dataset depression-filling algorithm. [Link](https://github.com/r-barnes/Barnes2016-ParallelPriorityFlood)
 * Large dataset flow accumulation algorithm. [Link](https://github.com/r-barnes/Barnes2016-ParallelFlowAccum)



Publications
============
The algorithms used in RichDEM have been published in the following articles. Every algorithm/program will provide its relevant citation information when run.

* Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion cell digital elevation models on desktops or clusters. Environmental Modelling & Software 92, 202–212. doi:[10.1016/j.envsoft.2017.02.022](https://doi.org/10.1016/j.envsoft.2017.02.022)

* Barnes, R., 2016. Parallel priority-flood depression filling for trillion cell digital elevation models on desktops or clusters. Computers & Geosciences 96, 56–68. doi:[10.1016/j.cageo.2016.07.001](https://doi.org/10.1016/j.cageo.2016.07.001)

* Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:[10.1016/j.cageo.2013.01.009](https://doi.org/10.1016/j.cageo.2013.01.009)

* Barnes, R., Lehman, C., Mulla, D., 2014b. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:[10.1016/j.cageo.2013.04.024](https://doi.org/10.1016/j.cageo.2013.04.024)

* Barnes, Lehman, Mulla. 2011. "Distributed Parallel D8 Up-Slope Area Calculation in Digital Elevation Models". Intn'l Conf. on Parallel & Distributed Processing Techniques & Applications. [Link](http://rbarnes.org/section/sci/2011_barnes_distributed.pdf)

* Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69, 14–47. doi:[10.1109/PROC.1981.11918](http://dx.doi.org/10.1109/PROC.1981.11918)

* Mulla et al. including Barnes. 2012. "Strategic Planning for Minnesota’s Natural and Artificial Watersheds". Report to the Minnesota LCCMR.

* Tarboton, D.G. 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Resources Research. Vol. 33. pp 309-319.

* Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. Earth surface processes and landforms 12, 47–56.

* Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of the Priority-Flood algorithm for filling depressions in raster digital elevation models. Computers & Geosciences 90, Part A, 87 – 96. doi:[10.1016/j.cageo.2016.02.021](http://dx.doi.org/10.1016/j.cageo.2016.02.021)



Credits
=======

RichDEM has been developed and tested using computational resources provided by
the [Minnesota Supercomputing Institute](https://www.msi.umn.edu/) (MSI) and the
U.S. National Science Foundation's [XSEDE](https://www.xsede.org/).

Funding for the development of RichDEM has been provided by the [Legislative-Citizen Commission on Minnesota Resources](http://www.lccmr.leg.mn/) (LCCMR), the U.S. National Science
Foundation [Graduate Research Fellowship](https://www.nsfgrfp.org/), and the U.S. Department of Energy
[Computational Science Graduate Fellowship](https://www.krellinst.org/csgf/).



Feedback
========

_If you see something, say something._

Users are encouraged to report any issues experienced with the code via Github's
issue tracker. Feedback is also accepted via email (rbarnes@umn.edu), though
this is highly discouraged as it does not provide a resource for others.
