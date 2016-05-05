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

    sudo apt-get install make openmpi-bin libgdal-dev libopenmpi-dev

If you wish (as I did) to compile the code on XSEDE, certain modules must be
loaded:

    module load intel/2015.2.164
    module load mvapich2_ib

Note that temporary files can be stored in:

    /oasis/scratch/comet/$USER/temp_project

or some similar directory.

Running `make` will produce an executable called `parallel_pf.exe`.



Running the Program
-------------------

`parallel_pf.exe` can be run without arguments from the command line to show a
comprehensive explanation of the program and its options. This same text is in
the file `help.txt`.

In order to process data, you will need to run `parallel_pf.exe` in MPI. For
example:

    mpirun -n 4 ./parallel_pf.exe one @offloadall dem.tif outroot -w 500 -h 500

In the foregoing example `-n 4` indicates that the program should be run in
parallel over four processes. One of these processes (the one with MPI rank #0)
acts as a master process. It does limited computation but stores information
from all of the other processes. This requires less memory than one would think,
as discussed in the manuscript.



Layout Files
------------

A layout file is a text file with the format:

    file1.tif, file2.tif, file3.tif,
    file4.tif, file5.tif, file6.tif, file7.tif
             , file8.tif,          ,

where each of fileX.tif is a tile of the larger DEM collectively described by
all of the files. All of fileX.tif must have the same shape; the layout file
specifies how fileX.tif are arranged in relation to each other in space.
Blanks between commas indicate that there is no tile there: the algorithm will
treat such gaps as places to route flow towards (as if they are oceans). Note
that the files need not have TIF format: they can be of any type which GDAL
can read. Paths to fileX.tif are taken to be relative to the layout file.

Several example layout files are included in the `tests/` directory and end with
the `.layout` extension.



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

    mpirun -output-filename timing -n 4 /usr/bin/time -v ./parallel_pf.exe one @offloadall dem.tif outroot -w 500 -h 500

This will store memory and timing information in files beginning with the stem
`timing`.



Testing
-------

For **running tests**, the following command will set you up on a Debian-based
system:

    sudo apt-get install python3-gdal python-gdal gdal-bin

The directory `tests` contains all of the information and layouts associated
with the tests described in the paper. The most immediately useful are probably
the `tests/beauford` test, which includes a small DEM suitable for testing the
correctness of various tile sizing configurations, and the `tests/srtm_small`
test (see the `README.md` file in that directory for further information), which
tests the "many" mode on a 3x3 excerpt of the SRTM Region 3 data.

Other subdirectories of `tests` are named for the dataset they pertain to and
contain directions for acquiring the datasets and example jobs for running them
using SLURM.

The `beauford` and `srtm_small` tests can be run using the `test.py` script.
This script can be running using one of the following: 

    ./test.py tests/beauford/beauford.tif
    ./test.py tests/srtm_small/srtm_small.layout

Once data has been acquired and placed in these directories.

In the case of a layout file being used, the `test.py` script will merge all of
the tiles together. This merged file, or, in the case of a single input file
being used, that file, will be depression filled using the algorithm in a
single-core mode. This generates an authoritative answer against which
correctness is checked. The program then iterates over many tile sizes to ensure
that they all compare correctly against this authoritative answer.



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