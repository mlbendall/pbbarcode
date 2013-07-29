.. pbbarcode Functional Specification
.. =======================================

.. Version


Introduction
````````````

This document describes the interface and input/output formats of the
``pbbarcode`` package command line tools. The package provides
utilities for annotating individual ZMWs directly from a bas.h5 file,
emitting fast[a|q] files for each barcode, labeling alignments stored
in a cmp.h5 file, and calling consensus on small amplicon sequencing
runs (requires pbdagcon)

At the moment, two analysis modes are supported. The ``symmetric``
mode supports barcode designs with two identical barcodes on both
sides of a SMRTbell, e.g., for barcodes (A, B), molecules are labeled
as A--A or B--B. The ``paired`` mode supports designs with two
distinct barcodes on each side of the molecule which are different,
but neither barcode appears without its mate. The minimum example is
given with the following barcodes: (ALeft, ARight, BLeft, BRight),
where the following barcode sets are checked: ALeft--ARight,
BLeft--BRight.

It is important to highlight that barcode fasta files specify a list
of available barcodes to evaluate. Depending on the scoring mode, the
barcodes are grouped together in different ways. For instance, in the
``symmetric`` case, the number of possible barcode outcomes are simply
the number of barcodes that are supplied to the routine in the fasta
file (see below for usage) plus an additional null barcode indicating
that no barcode could be evaluated (denoted by: '--'). Labels like
this (A--A) are used in the final outputs. 

The software is implemented as a standard python package. Barcodes are
labeled according to the following high-level logic. For each
molecule, all adapters are found. For each adapter, we align (using
standard Smith-Watterman alignment) each barcode and its reverse
complement to flanking sequence of the adapter. We record the score
for this barcode. Depending on the ``mode``, we then determine which
barcode(s) are maximally scoring. We store the two maximally scoring
barcodes, the sum of their alignment scores across the adapters. The
average barcode score then can be given approximately by:
score/(2*number-of-adapters); each adapter yields 2 barcode
alignments.


Input and output
````````````````

labelZMWs
---------

| usage: pbbarcode.py labelZMWs 
|                                [--outDir OUTDIR]   
|                                [--outFofn OUTFOFN]  
|                                [--adapterSidePad ADAPTERSIDEPAD] 
|                                [--insertSidePad INSERTSIDEPAD] 
|                                [--scoreMode {symmetric,paired}]
|            barcode.fasta input.fofn

The ``labelZMWs`` command takes an input.fofn representing a set of
bas.h5 files to operate on. Additionally, it takes a barcode.fasta
file. Depending on ``scoreMode``, the fasta file will be processed in
different ways. Specifically, in ``paired`` mode, each two consecutive
barcodes in the file are considered a set.

The parameters, ``adapterSidePad`` and ``insertSidePad`` represents
how many bases should be considered on each side of the putative
barcode. These parameters are constrained such that:
``|adapterSidePad| + |insertSidePad| + |barcode| < 64``.

Finally, users have the option to specify a different output location
for the various outputs. Specifically, for each bas.h5 file in
input.fofn, a bc.h5 (barcode hdf5) file is generated. These files are
listed in the file ``outFofn`` which is typically just called
``barcode.fofn``. See below for a description of the barcode hdf5
file.

labelAlignments
---------------

| usage: pbbarcode.py labelAlignments barcode.fofn aligned_reads.cmp.h5

The ``labelAlignments`` command takes as input a barcode.fofn computed
from a call to ``labelZMWs`` and a cmp.h5 file where the barcode
information is written to. See below for a description of the cmp.h5
file additions. 


emitFastqs
----------
| usage: pbbarcode.py emitFastqs [--outDir output.dir] 
|                                [--trim TRIM]
| 			         [--subreads]
|                         input.fofn barcode.fofn

The ``emitFastqs`` command takes as input both an input.fofn for the
bas.h5 files as well as a barcode.fofn from a call to labelZMWs. The
optional parameter ``outDir`` dictates where the files will be written
to. For each detected barcode, a fast[a|q] file will be emitted with
all of the reads for that barcode. The ``trim`` parameter dictates how
much of the read should be trimmed off. The default parameter for
``trim`` is the length of the barcode (which is stored in the barcode
hdf5 files). At the moment, all barcodes in the barcode fasta file
must be the same length, therefore only a constant trim value is
supported. In practice, one can aggressively trim in order to ensure
that extra bases aren't left on the ends of reads. Finally, the
``subreads`` parameter dictates whether subreads or CCS reads should
be returned with the default being CCS reads. The ``subreads``
parameter is only inspected if the input.fofn contains both CCS and
subread data, if the input.fofn contains only subread or CCS data then
that is returned irrespective of the state of the the ``subreads``
parameter.

Dependencies
````````````
The pbbarcode package depends on a standard pbcore installation.

Barcode HDF5 File
`````````````````

The barcode hdf5 file, ``bc.h5``, represents a simple data store for
barcode calls and their scores for each ZMW. Generally, a user need
not interact with barcode hdf5 file, but can use the results stored in
either the resulting cmp.h5 file or fast[a|q] files. The barcode hdf5
file contains the following structure:

/BarcodeCalls/best - (nZMWs, 6)[32-bit integer] dataset with the
following columns: 

``holeNumber,nAdapters,barcodeIdx1,barcodeScore1,barcodeIdx2,barcodeScore2``

Additionally, the ``best`` dataset has the following attributes:

|  movieName   = m120408_042614_richard_c100309392550000001523011508061222_s1_p0
|  columnNames = holeNumber,nAdapters,barcodeIdx1,barcodeScore1,barcodeIdx2,barcodeScore2
|  scoreMode   = symmetric
|  barcodes    = bc_1, ... 

The two barcodeIdx1 and barcodeIdx2 columns are indices into
``barcodes`` attribute. The ``scoreMode`` is scoring mode used to
align the barcodes. The ``barcodes`` attribute correspond to the
barcode.fasta file names. 

Additions to the compare HDF5 (cmp.h5) File
```````````````````````````````````````````

In addition to the barcode hdf5 file, a call to ``labelAlignments``
will annotate a cmp.h5 file. This annotation is stored in ways
consistent with the cmp.h5 file format. Specifically, a new group: 

| /BarcodeInfo/
|   ID   (nBarcodeLabels + 1, 1)[32-bit integer] 
|   Name (nBarcodeLabels + 1, 1)[variable length string]

In addition to the /BarcodeInfo/ group, the key dataset which assigns
alignments to barcodes is located at:

/AlnInfo/Barcode (nAlignments, 3)[32-bit integer] with the following
colums:

``index,count,bestIndex,bestScore,secondBestIndex,secondBestScore``

Here index refers to the index into the ``Name`` vector, score
corresponds to the sum of the scores for the barcodes, and finally,
count refers to the number of adapters found in the molecule.
