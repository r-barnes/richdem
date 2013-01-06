RichDem
=======

Author: Richard Barnes (rbarnes@umn.edu)

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly process even very large DEMs.

RichDEM can use both the D8 and D-infinite (Tarboton) flow metrics. It can resolve terrain depressions (or pits) either by filling or by channel carving. It can calculate contributing/up-slope areas, slopes, curvatures, and aspects.

Future versions may use Intel's Thread Building Blocks to achieve additional increases in speed.

To use RichDEM, edit "main.cpp" to suit your needs and then simply run "make".

To generate DOxygen documentation use "doxygen doxy.conf".

Requires:
OpenMP
