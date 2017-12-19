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
