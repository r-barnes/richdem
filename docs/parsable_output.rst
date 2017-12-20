Parsable Output
===================

Every line of output from RichDEM begins with one of the following characters,
making it easy to parse with a machine.

 ======== ===================================================================
 Tag      Meaning
 ======== ===================================================================
 **A**    Algorithm name
 **a**    Analysis command: the command line used to run the program
 **c**    Configuration information: program version, input files, and command line options, &c.
 **C**    Citation for algorithm
 **d**    Debugging info
 **E**    Indicates an error condition
 **i**    I/O: Amount of data loaded from disk
 **m**    Miscallaneous counts
 **n**    I/O: Amount of data transferred through a network
 **p**    Progress information: inform the user to keep calm because we're carrying on.
 **r**    Amount of RAM used
 **t**    Timing information: How long stuff took
 **W**    Indicates a warning
 ======== ===================================================================

All output data shall have the form:

    <INDICATOR TAG> <MESSAGE/MEASUREMENT NAME> [= <VALUE> [UNIT]]

The amount of whitespace may very for aesthetic purposes.
