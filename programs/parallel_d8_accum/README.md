Barnes2016-ParallelFlowAccum
============================

**Title of Manuscript**:
Parallel Non-divergent Flow Accumulation For Trillion Cell Digital Elevation Models On Desktops Or Clusters

**Authors**: Richard Barnes

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**DOI Number of Manuscript**:
[10.1016/j.envsoft.2017.02.022](http://dx.doi.org/10.1016/j.envsoft.2017.02.022)

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2016-ParallelFlowAccum)

This repository contains a reference implementation of the algorithms presented
in the manuscript above, along with information on acquiring the various
datasets used, and code to perform correctness tests.



Abstract
--------

Continent-scale datasets challenge hydrological algorithms for processing
digital elevation models. Flow accumulation is an important input for many such
algorithms; here, I parallelize its calculation. The new algorithm works on one
or many cores, or multiple machines, and can take advantage of large memories or
cope with small ones. Unlike previous algorithms, the new algorithm guarantees a
fixed number of memory access and communication events per raster cell. In
testing, the new algorithm ran faster and used fewer resources than previous
algorithms exhibiting ∼30% strong and weak scaling efficiencies up to 48 cores
and linear scaling across datasets ranging over three orders of magnitude. The
largest dataset tested has two trillion (2 · 10^12) cells. With 48 cores,
processing required 24 minutes wall-time (14.5 compute-hours). This test is
three orders of magnitude larger than any previously performed in the
literature. Complete, well-commented source code and correctness tests are
available for download from Github.



Compilation
-----------

The program is compiled as part of the RichDEM. Please look at RichDEM's main
README.md file for compilation instructions.

You must also obtain certain prerequisites:

    sudo apt install make openmpi-bin libgdal-dev libopenmpi-dev

If you wish (as I did) to compile the code on XSEDE, certain modules must be
loaded:

    module load intel/2015.2.164
    module load mvapich2_ib

To compile the programs run:

    make

The result is a program called `parallel_d8_accum.exe`.

Running the above compiles the program to run the _cache_ strategy. Using `make
compile_with_compression` will enable the _cacheC_ strategy instead. This
strategy is not compiled by default because it requires the Boost Iostreams
library. This libary can be installed with:

    sudo apt-get install libboost-iostreams-dev



Further details
---------------

For further details on testing, layout files, and so on, please see the source code's
[README.md](https://github.com/r-barnes/richdem/blob/master/programs/parallel_d8_accum/README.md).



MPI Profiling
-------------
Although the program tracks its total communication load internally, I have also
used [mpiP](http://mpip.sourceforge.net/) to profile the code's communication.
The code can be downloaded [here](http://mpip.sourceforge.net/) and compiled
with:

    ./configure --with-binutils-dir=/usr/lib
    make shared
    make install #Installs to a subdirectory of mpiP

Prerequisites include: `binutils-dev`.

mpiP can be used to profile _any_ MPI program without the need to compile it
with the program. To do so, run the following line immediately before launching
`mpirun`:

    export LD_PRELOAD=path/to/libmpiP.so

Although the program tracks its maximum memory requirements internally, I have
also used `/usr/bin/time` to record this. An example of such an invocation is:

    mpirun -output-filename timing -n 4 /usr/bin/time -v ./parallel_d8_accum.exe one @offloadall dem.tif outroot -w 500 -h 500

This will store memory and timing information in files beginning with the stem
`timing`.



RichDEM
-------

This code is part of the RichDEM codebase, which includes state of the art
algorithms for quickly performing hydrologic calculations on raster digital
elevation models. The full codebase is available at
[https://github.com/r-barnes](https://github.com/r-barnes)
