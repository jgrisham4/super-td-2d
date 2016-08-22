super-td-2d
=========================
super-td-2d stands for supersonic triple-deck two-dimensional.  The code is an
implementation of a triple-deck solver described by Cassel, et al., in their
JFM paper from 1996.  The code solves the nonlinear triple-deck problem using
second-order accurate finite differences in space and a first-order accurate 
temporal discretization.  Results are written out in ASCII Tecplot fepoint 
format.  Results such as scaled shear stress and scaled pressure at the wall
are also written out to an ASCII text file.  Inputs are given through a 
token-based input file which is parsed at run time.  Restart capabilities are
also built-in using the Boost serialization library.

Prerequisites:
------------------------
The only dependency of this code is Boost.  The Boost serialization library
is used to write objects to files which are then used for checkpointing and
restarting.  I have tested this code with Boost 1.60.0 on my Arch Linux machine
and Boost 1.56.0 on my Fedora machine.  Both work fine.

Installation:
------------------------
The code can be built by editing the makefile in this directory.  Make sure
the linker can find the libboost_serialization.so or .a library.  The executable
will be built in src.  

License:
-----------------------
Distributed under GNU GPLv3.
