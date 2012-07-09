
.. pbbarcode Functional Specification
.. =======================================

.. Version


Introduction
````````````
This document describes the interface and input/output formats of the
``pbbarcode`` package command line tools. The package provides
utilities for annotating ZMWs as well as alignments stored in
``bas.h5`` and ``cmp.h5`` files respectively. At the moment, two modes
are supported. The ``symmetric`` mode supports barcode designs with
two identical barcodes on both sides of a SMRTbell. The ``asymmetric``
mode supports barcode designs with different barcodes on each side.

Software Overview
-----------------


Functional Requirements
-----------------------



Input and output
````````````````

``pbbarcode.py labelMovie`` requires the following input files:
input
    - input.bas.h5: Input base hdf5 file.
    - barcode.fasta : FASTA file with barcode sequences.
output
    - barcode.h5 file : The barcode h5 file is described in section ***

``pbbarcode.py labelAlignments`` requires the following input files:
    - barcode.fofn

.. note::
        **Input cmp.h5 file requirements**



Dependencies
````````````
The pbbarcode package depends on a pbcore installation.

Barcode HDF5 File
`````````````````

The barcode hdf5 file, ``bc.h5``, represents a simple data store for
barcode calls and scores. 

