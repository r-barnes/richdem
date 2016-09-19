Hierarchy of Correctness
========================

RichDEM uses a few different methodologies to establish the correctness of its
algorithms and implementations. In the simplest cases, correctness is
established via inspection of the code. In the most complicated cases, the
output of a simpler algorithm whose correctness can be more easily established
is compared against the output of a more complex algorithm.

The following describes how the correctness of each algorithm is established.

* `d8_flow_accum()`: Correctness is established via inspection and hand-made
                     tests in `tests/flow_accum`.

* `rd_flow_accumulation.exe`: Correctness is established because this program
                              wraps `d8_flow_accum()`.

* `parallel_d8_accum.exe`: Correctness is established via comparison against
                           hand-made tests in `tests/flow_accum` and by 
                           comparing outputs on larger/realistic inputs against
                           `rd_flow_accumulation.exe`.




Tests (TODO)
============

These directories contain information useful in acquiring and performing the
tests described in the manuscript. See `README.md` files in individual
directories for more information.

Several utility programs are included:

 * `find_square.py`: Finds the largest contiguous square in a layout file and 
                     prints appropriate layout files

 * `display_layout.py`: Displays a layout file to the terminal for verification
