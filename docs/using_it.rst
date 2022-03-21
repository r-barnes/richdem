Ways To Use It
==============

Python package from PyPI
------------------------

Get the package with:

.. code-block:: bash

    pip3 install richdem

And use:

.. code-block:: python

    import richdem

The command:

    help(richdem)

provides all the relevant documentation.

Python package from source
--------------------------

Enter the `wrappers/pyrichdem` directory and run:

.. code-block:: bash

    python3 setup.py install --user



As A Command-line Tool
----------------------

To get the command-line tools, install the Python package with:

.. code-block:: bash

    pip3 install richdem

The command-line tools are all named `rd_*`, so typing `rd_` on your command-
line and hitting tab a few times should give you the full list of what's
available.



As A Library
------------

As an overview, upon compilation, point your library search path to the
`richdem/include` directory. Include various files using, e.g.

.. code-block:: c++

    #include "richdem/common/Array2D.hpp"

All files include extensive documentation. At this stage the location of certain
functions may be subject to change. This will be noted in the `NEWS` file. (TODO)

More concretely, there are a number of compilation options to consider. With the
GCC compiler the flag `FLAG` is activated by passing the `-DFLAG` command-line
argument.

 * `NOPROGRESS` turns off progress bars

 * `RICHDEM_DEBUG` turns on line numbers and filenames in RichDEM's output

 * `RICHDEM_LOGGING` turns on outputs such as notices about memory allocation,
   intermediate products, various progress indicators, and so on. Enabling this
   requires inclusion of the `richdem.cpp` file.

 * `RICHDEM_GIT_HASH`. this should be set to the git hash of the code you checked
   out. A value of `RICHDEM_GIT_HASH=$(git rev-parse HEAD)` is usually good.

 * `RICHDEM_COMPILE_TIME`. Date and time the code was compiled. A value of
   `RICHDEM_COMPILE_TIME=$(date -u +'%Y-%m-%d %H:%M:%S UTC')` is usually good.

 * `USEGDAL`. Indicates that GDAL functionality should be included in the
   library. This allows reading/writing rasters from various file types. It also
   complicates compilation slightly, as discussed below.

 * `NDEBUG` turns off a bunch of range-checking stuff included in the standard
   library. Increases speed slightly, butm akes debugging crashes and such more
   difficult.

Setting up compilation works like this:

.. code-block:: bash

    CXXFLAGS="--std=c++17 -g -O3 -Wall -Wno-unknown-pragmas -Irichdem/include"
    CXXFLAGS="$CXXFLAGS -DRICHDEM_LOGGING"

C++17 or higher is necessary to compile. Include other RichDEM flags as desired.
Note that the `-g` flag doesn't slow things down, though it does increase the
size of your executable somewhat. It's inclusion is always recommended for
anything other than distributed production code because it makes debugging much
easier. The `-O3` flag should be replaced by an optimization level or set of
your choice. `-Wno-unknown-pragmas` hides warning messages from OpenMP if you
choose not to compile with it. `-Wall` produces many helpful warning messages;
compiling without `-Wall` is foolish. `-Irichdem/include` connects your code
with RichDEM.

If you plan to use GDAL, include the following:

.. code-block:: bash

    GDAL_LIBS="`gdal-config --libs`"
    GDAL_CFLAGS="`gdal-config --cflags` -DUSEGDAL"
    LIBS="$GDAL_LIBS"
    CXXFLAGS="$CXXFLAGS $GDAL_CFLAGS"

If you plan to use RichDEM's parallel features include the following:

    LIBS="$LIBS -fopenmp"

Finally, put it all together:

.. code-block:: bash

    g++ $CXXFLAGS -o my_program.exe my_program.cpp $LIBS



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
