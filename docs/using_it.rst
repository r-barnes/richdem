
Using It
========

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

