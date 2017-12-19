Design Philosophy
=================

The design of RichDEM is guided by these principles:

- **Algorithms will be well-tested.** Every algorithm is verified by a rigorous
  testing procedure. See below.

- **Algorithms will be fast, without compromising safety and accuracy.** The
  algorithms used in RichDEM are state of the art, permitting analyses that
  would take days on other systems to be performed in hours, or even minutes.

- **Algorithms will be available as libraries, whenever possible.** RichDEM is
  designed as a set of header-only C++ libraries, making it easy to include in
  your projects and easy to incorporate into other programming languages.
  RichDEM also includes apps, which are simple wrappers around the algorithms, 
  and a limited, but growing, set of algorithms which may have special
  requirements, like MPI, that make them unsuitable as libraries. These are 
  available as programs.

- **Programs will have a command-line interface, not a GUI.** Command-line
  interfaces are simple to use and offer extreme flexibility for both users and
  programmers. They are available on every type of operating system. RichDEM
  does not officially support any GUI. Per the above, encapsulating RichDEM in
  a high-level interface of your own is not difficult.

- **Algorithms and programs will be portable.** Linux, Mac, and Windows should
  all be supported.

- **The code will be beautiful.** RichDEM's code utilizes sensible variable
  names and reasonable abstractions to make it easy to understand, use, and
  design algorithms. The code contains extensive internal documentation which is
  DOxygen compatible.

- **Programs and algorithms will provide useful feedback.** Progress bars will 
  appear if desired and the output will be optimized for machine parsing.

- **Analyses will be reproducible.** Every time you run a RichDEM command that
  command is logged and timestamped in the output data, along with the version
  of the program you created the output with. Additionally, a history of all
  previous manipulations to the data is kept. Use `rd_view_processing_history`
  to see this.**
