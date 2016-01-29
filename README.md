Barnes2016-DistributedPriorityFlood
===================================

**Title of Manuscript**:
Flexibly Distributed Priority-Flood Depression Filling For Teracell DEMs

**Authors**: Richard Barnes

**Corresponding Author**: Richard Barnes (rbarnes@umn.edu)

**DOI Number of Manuscript**
TODO

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2016-DistributedPriorityFlood)
 * [Journal's GitHub Repository](TODO)

This repository contains a reference implementation of the algorithms presented
in the manuscript above, along with information on acquiring the various
datasets used, and code to perform correctness tests.




Abstract
--------
High-resolution digital elevation models (DEMs) are increasingly available;
however, attempts to use these in existing hydroanalysis algorithms has created
problems which take too long to solve and which are too big to fit into a
computer's working memory. This has led to research on parallel algorithms and
algorithms which explicitly manage memory. Parallel approaches do not scale well
because they require many nodes and frequent communication between these nodes.
Memory-managing algorithms have to read and write subdivisions of the data
numerous times because they suffer from low access locality. Here, I adopt a
tile-based approach which unifies the parallel and memory-managing paradigms. I
show that it is possible to **perform depression-filling on a tiled DEM with a
fixed number of memory access and communication events per tile**, regardless of
the size of the DEM. The result is an algorithm that works equally well on one
core, multiple cores, or multiple machines and can take advantage of large
memories or cope with small ones. The largest dataset on which I run the
algorithm has 2 trillion (2*10^12) cells. With 48 cores, processing required 291
minutes to completion and 9.4 compute-days. This test is three orders of
magnitude larger than any previously performed in the literature, but took only
1-2 orders of magnitude more compute-time. Complete, well-commented source code
and correctness tests are available for download from a repository.





Compilation
-----------

The program can be produced simply by running **make**. However, certain
prerequisites are necessary for this to be successful.

For **compilation**, the following command will set you up on a Debian-based
system:

    sudo apt-get install make openmpi-bin libgdal-dev libopenmpi-dev \
                         libboost-mpi-dev libboost-filesystem-dev    \
                         libboost-system-dev 

If you wish (as I did) to compile the code on XSEDE, certain modules must be
loaded:

    module unload intel/2013_sp1.2.144
    module load boost/1.55.0
    module load intel/2015.2.164
    module load mvapich2_ib

Note that temporary files can be stored in:

    /oasis/scratch/comet/$USER/temp_project

or some similar directory.




Running the Algorithm
---------------------

Running `make` will produce an executable called `parallel_pf.exe`. This can be
run without arguments from the command line to show a comprehensive explanation
of the program and its options. This same text is in the file `help.txt`.




Testing
-------

For **running tests**, the following command will set you up on a Debian-based
system:

    sudo apt-get install python3-gdal python-gdal gdal-bin

The directory `tests` contains all of the information and layouts associated
with the tests described in the paper. The most immediately useful are probably
the `tests/beauford` test, which includes a small DEM suitable for testing the
correctness of various tile sizing configurations, and the
`tests/srtm_region/nasa_srtm3_small.layout` test (see the `README.md` file in
that directory for further information), which tests the "many" mode on a 3x3
excerpt of the SRTM Region 3 data.





RichDEM
-------

This code is part of the RichDEM codebase, which includes state of the art
algorithms for quickly performing hydrologic calculations on raster digital
elevation models. The full codebase is available at
[https://github.com/r-barnes](https://github.com/r-barnes)




TODO
----

Different MPI Polling methods
https://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu