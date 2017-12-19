RichDEM
=======

Author: Richard Barnes (rbarnes@umn.edu)

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly
process even very large DEMs.

RichDEM offers a variety of flow metrics, such as D8 and D∞. It can flood or
breach depressions. It can calculate flow accumulation, slops, curvatures, &c.

Please cite RichDEM (see below).



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

* **Analyses will be reproducible." Every time you run a RichDEM command that
  command is logged and timestamped in the output data, along with the version
  of the program you created the output with. Additionally, a history of all
  previous manipulations to the data is kept. Use `rd_view_processing_history`
  to see this.



Using It
========



Citing It
---------

As of 883ea734e957, David A. Wheeler's SLOCCount estimates the value of RichDEM
at $240,481 and 1.78 person-years of development effort. This value is yours to
use, but citations are encouraged.

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

Users are also encouraged to cite the particular algorithms used. Citations to
these will be printed whenever an app, program, or library function is run.
Although I have written all of the code in this library, some of the algorithms
were discovered or invented by others, and they deserve credit for their good
work.



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



Testing Methodology
===================

TODO



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
          carrying on.

 * **m**: Miscallaneous counts

 * **n**: I/O: Amount of data transferred through a network

 * **p**: Progress information: inform the user to keep calm because we're
 
 * **r**: Amount of RAM used

 * **t**: Timing information: How long stuff took

 * **W**: Indicates a warning


All output data shall have the form:

    <INDICATOR CHARACTER> <MESSAGE/MEASUREMENT NAME> [= <VALUE> [UNIT]]

The amount of whitespace may very for aesthetic purposes.






Publications
============
The algorithms used in RichDEM have been published in the following articles. Every algorithm/program will provide its relevant citation information when run.

* Barnes, R., 2016. Non-divergent flow accumulation for trillion cell digital elevation models on desktops or clusters. In Review.

* Barnes, R., 2016. Parallel priority-flood depression filling for trillion cell digital elevation models on desktops or clusters. Computers & Geosciences. doi:[10.1016/j.cageo.2016.07.001](http://dx.doi.org/10.1016/j.cageo.2016.07.001)

* Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of the Priority-Flood algorithm for filling depressions in raster digital elevation models. Computers & Geosciences 90, Part A, 87 – 96. doi:http://dx.doi.org/10.1016/j.cageo.2016.02.021

* Barnes, Lehman, Mulla. 2013. "An Efficient Assignment of Drainage Direction Over Flat Surfaces in Raster Digital Elevation Models". Computers &amp; Geosciences. doi: [10.1016/j.cageo.2013.01.009](http://dx.doi.org/10.1016/j.cageo.2013.01.009)

* Barnes, Lehman, Mulla. 2013. "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models". Computers &amp; Geosciences. doi: [10.1016/j.cageo.2013.04.024](http://dx.doi.org/10.1016/j.cageo.2013.04.024)

* Mulla et al. including Barnes. 2012. "Strategic Planning for Minnesota’s Natural and Artificial Watersheds". Report to the Minnesota LCCMR.

* Barnes, Lehman, Mulla. 2011. "Distributed Parallel D8 Up-Slope Area Calculation in Digital Elevation Models". Intn'l Conf. on Parallel & Distributed Processing Techniques & Applications. [Link](http://rbarnes.org/section/sci/2011_barnes_distributed.pdf)

* Tarboton, D.G. 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Resources Research. Vol. 33. pp 309-319.

* Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69, 14–47. doi:10.1109/PROC.1981.11918

* Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. Earth surface processes and landforms 12, 47–56.



Credits
=======

RichDEM has been developed and tested using computational resources provided by
the [Minnesota Supercomputing Institute][1] (MSI) and the U.S. National Science
Foundation's [XSEDE][2].

Funding for the development of RichDEM has been provided by the [Legislative-Citizen Commission on Minnesota Resources][5] (LCCMR), the U.S. National Science
Foundation [Graduate Research Fellowship][3], and the U.S. Department of Energy
[Computational Science Graduate Fellowship][4].



Feedback
========

_If you see something, say something._

Users are encouraged to report any issues experienced with the code via Github's
issue tracker. Feedback is also accepted via email (rbarnes@umn.edu), though
this should be used only in cases wherein the issue tracker is not available.




[1]: https://www.msi.umn.edu/
[2]: https://www.xsede.org/
[3]: https://www.nsfgrfp.org/
[4]: https://www.krellinst.org/csgf/
[5]: http://www.lccmr.leg.mn/