2022-02-11 (2.3.1)
==================

Fix pyrichdem to import all cpp source files correctly
Clean up pyrichdem setup.py
Add .log files to .gitignore

2022-02-09 (2.3.0)
==================

Drop PyBind11 in favour of using pypi's version
Update to doctest v2.4.8

2018-07-13 (2.2.9)
==================

Changed `isnan()` logic in `grid_cell.hpp` to facilitate MSVC compilation.


2018-06-21 (2.2.8)
==================

Release for Zenodo DOI.


2018-06-03 (2.2.7)
==================

Fixes some compilation issues for Windows caused by MSVC not having M_PI by default.


2018-06-02 (2.2.6)
==================

Fixes some compilation issues for Windows.


2018-05-31 (2.2.5)
==================

Travis-CI deployment to Github works.


2018-05-31 (2.2.0)
==================

Includes D4 flow routing, depression filling, and depression breaching
Uses `chrono` for timing instead of `sys/time.h`
Includes Windows compatibility improvements
