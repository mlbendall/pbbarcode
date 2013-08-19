Overview of the pbbarcode package
=================================

The *pbbarcode* package provides tools for annotating PacBio
sequencing reads with barcode information. Typically, *pbbarcode.py*
is called in context of a SMRTPipe workflow as opposed to directly on
the command line, however, users are encouraged to utilize the
command-line utility directly, as more options are available.  

The *pbbarcode* package provides a multi-command line tool
*pbbarcode.py* which currently has the following sub-commands:  

* labelZmws
* labelAlignments
* emitFastqs
* consensus

The first three sub-commands depend on only *pbcore* and its
dependencies, the fourth, *consensus*, depends on the *pbdagcon*
package and is considered experimental.  

For more details on the package, please see docs/index.rst for more
information.

Installation
============

Typically, the *pbbarcode* package is installed within an installation
of SMRTPipe, however, it can be installed by itself using::

   make install

To test that everything is installed correctly, one should
additionally issue a::

   make test
