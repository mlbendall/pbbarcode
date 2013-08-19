.. pbbarcode Functional Specification
.. =======================================

.. Version


Introduction
````````````
This document describes the interface and input/output formats of the
``pbbarcode`` package command line tools. The package provides
utilities for annotating individual ZMWs directly from a bas.h5 file,
emitting fast[a|q] files for each barcode, labeling alignments stored
in a cmp.h5 file, and calling consensus on small amplicons (requires
``pbdagcon``)

At the moment, Barcodes can be scored in two different ways:
``symmetric`` and ``paired``. Symmetric mode supports barcode designs
with two identical barcodes on both sides of a SMRTbell, e.g., for
barcodes (A, B), molecules are labeled as A--A or B--B. The ``paired``
mode supports designs with two distinct barcodes on each side of the
molecule, but neither barcode appears without its mate. The minimum
example is given with the following barcodes: (ALeft, ARight, BLeft,
BRight), where the following barcode sets are checked: ALeft--ARight,
BLeft--BRight.

It is important to highlight that a barcode FASTA file specifies a
list of available barcodes to evaluate. Depending on the scoring mode,
the barcodes are grouped together in different ways. For instance, in
the ``symmetric`` case, the number of possible barcode outcomes are
simply the number of barcodes that are supplied to the routine in the
FASTA file (see below for usage) plus an additional ``NULL`` barcode
indicating that no barcode could be evaluated (denoted by:
'--'). Labels like this (A--A) are used in the final outputs. In the
``paired`` mode, the number of possible barcode outcomes are half the
number of the sequences in the FASTA file plus the ``NULL``
barcode. The ``NULL`` barcode indicates that no attempt was made to
score the molecule or it was filtered out by the user's criteria. The
majority of cases when a molecule is not scored are related to not
observing any adapters. If a user has executed a "hot-start" run, the
user can try the '--scoreFirst' parameter to attempt to label the
first adapter's barcode. This increases the yield of the labeleing
procedure at the expense of some probably false positives. 

The software is implemented as a standard python package. Barcodes are
labeled according to the following high-level logic. For each
molecule, all adapters are found. For each adapter, we align (using
standard Smith-Watterman alignment) each barcode and its reverse
complement to flanking sequence of the adapter. If two complete
flanking sequences are available, we divide by 2, else 1 if only one
flanking sequence was available (average score at adapter). This
allows the scores across adapters to be on the same scale (chimera
detection). Depending on the ``mode``, we then determine which
barcode(s) are maximally scoring. We store the two maximally scoring
barcodes, the sum of their alignment scores across the adapters. The
average barcode score then can be given approximately by:
total-score/number-of-adapters. At the moment, the alignment
parameters are fixed at:


.. table:: SW Match Parameters
+----------+----------+
|type      |score     |
|          |          |
+----------+----------+
|insertion |-1        |
|          |          |
+----------+----------+
|deletion  |-1        |
|          |          |
+----------+----------+
|missmatch |-2        |
|          |          |
+----------+----------+
|match     |2         |
|          |          |
+----------+----------+

Input and output
````````````````

labelZmws
---------
  usage: pbbarcode.py labelZmws [-h] [--outDir OUTDIR] [--outFofn OUTFOFN]
                                [--adapterSidePad ADAPTERSIDEPAD]
                                [--insertSidePad INSERTSIDEPAD]
                                [--scoreMode {symmetric,paired}]
                                [--maxAdapters MAXADAPTERS] [--scoreFirst]
                                [--startTimeCutoff STARTTIMECUTOFF]
                                [--nZmws NZMWS] [--nProcs NPROCS]
                                [--saveExtendedInfo]
                                barcode.fasta input.fofn

  Creates a barcode.h5 file from base h5 files.

  positional arguments:
    barcode.fasta         Input barcode fasta file
    input.fofn            Input base fofn

  optional arguments:
    -h, --help            show this help message and exit
    --outDir OUTDIR       Where to write the newly created barcode.h5 files.
                          (default: /home/UNIXHOME/jbullard/projects/software/bi
                         oinformatics/tools/pbbarcode/doc)
    --outFofn OUTFOFN     Write to outFofn (default: barcode.fofn)
    --adapterSidePad ADAPTERSIDEPAD
                          Pad with adapterSidePad bases (default: 4)
    --insertSidePad INSERTSIDEPAD
                          Pad with insertSidePad bases (default: 4)
    --scoreMode {symmetric,paired}
                          The mode in which the barcodes should be scored.
                          (default: symmetric)
    --maxAdapters MAXADAPTERS
                          Only score the first maxAdapters (default: 20)
    --scoreFirst          Whether to try to score the leftmost barcode in a
                          trace. (default: False)
    --startTimeCutoff STARTTIMECUTOFF
                          Reads must start before this value in order to be
                          included when scoreFirst is set. (default: 10.0)
    --nZmws NZMWS         Use the first n ZMWs for testing (default: -1)
    --nProcs NPROCS       How many processes to use (default: 8)
    --saveExtendedInfo    Whether to save extended information tothe barcode.h5
                          files; this information is useful for debugging and
                                                  chimera detection (default: False)

The ``labelZmws`` command takes an input.fofn representing a set of
bas.h5 files to operate on. Additionally, it takes a barcode.fasta
file. Depending on ``scoreMode``, the FASTA file will be processed in
different ways. Specifically, in ``paired`` mode, each two consecutive
barcodes in the file are considered a set.

The parameters, ``adapterSidePad`` and ``insertSidePad`` represents
how many bases should be considered on each side of the putative
barcode. These parameters are constrained such that:
``|adapterSidePad| + |insertSidePad| + |barcode| < 65``.

Users have the option to specify a different output location
for the various outputs. Specifically, for each bas.h5 file in
input.fofn, a bc.h5 (barcode hdf5) file is generated. These files are
listed in the file ``outFofn`` which is typically just called
``barcode.fofn``. See below for a description of the barcode hdf5
file.


labelAlignments
---------------
  usage: pbbarcode.py labelAlignments [-h]
                                      [--minAvgBarcodeScore MINAVGBARCODESCORE]
                                      [--minNumBarcodes MINNUMBARCODES]
                                      [--minScoreRatio MINSCORERATIO]
                                      barcode.fofn aligned_reads.cmp.h5

  Adds information about barcode alignments to a cmp.h5 file from a previous
  call to "labelZmws".

  positional arguments:
    barcode.fofn          input barcode fofn file
    aligned_reads.cmp.h5  cmp.h5 file to add barcode labels

  optional arguments:
    -h, --help            show this help message and exit
    --minAvgBarcodeScore MINAVGBARCODESCORE
                          ZMW Filter: exclude ZMW if average barcode score is
                          less than this value (default: 0.0)
    --minNumBarcodes MINNUMBARCODES
                          ZMW Filter: exclude ZMW if number of barcodes observed
                          is less than this value (default: 1)
    --minScoreRatio MINSCORERATIO
                          ZMW Filter: exclude ZMWs whose best score divided by
                          the 2nd best score is less than this ratio (default:
                          1.0)
                          

The ``labelAlignments`` command takes as input a barcode.fofn computed
from a call to ``labelZMWs`` and a cmp.h5 file where the barcode
information is written to. See below for a description of the cmp.h5
file additions.  



emitFastqs
----------
  usage: pbbarcode.py emitFastqs [-h] [--outDir output.dir] [--subreads]
                                 [--unlabeledZmws] [--trim TRIM] [--fasta]
                                 [--minMaxInsertLength MINMAXINSERTLENGTH]
                                 [--hqStartTime HQSTARTTIME]
                                 [--minReadScore MINREADSCORE]
                                 [--minAvgBarcodeScore MINAVGBARCODESCORE]
                                 [--minNumBarcodes MINNUMBARCODES]
                                 [--minScoreRatio MINSCORERATIO]
                                 input.fofn barcode.fofn

  Takes a bas.h5 fofn and a barcode.h5 fofn and produces a fast[a|q] file for
  each barcode.

  positional arguments:
    input.fofn            input base or CCS fofn file
    barcode.fofn          input barcode.h5 fofn file

  optional arguments:
    -h, --help            show this help message and exit
    --outDir output.dir   output directory to write fastq files (default: /home/
                          UNIXHOME/jbullard/projects/software/bioinformatics/too
                          ls/pbbarcode/doc)
    --subreads            whether to produce fastq files for the subreads;the
                          default is to use the CCS reads. This option
                          onlyapplies when input.fofn has both consensus and raw
                          reads,otherwise the read type from input.fofn will be
                          returned. (default: False)
    --unlabeledZmws       whether to emit a fastq file for the unlabeled ZMWs.
                          These are the ZMWs where no adapters are found
                          typically (default: False)
    --trim TRIM           trim off barcodes and any excess constant sequence
                          (default: 20)
    --fasta               whether the files produced should be FASTA files
                          asopposed to FASTQ (default: False)
    --minMaxInsertLength MINMAXINSERTLENGTH
                          ZMW Filter: exclude ZMW if the longest subreadis less
                          than this amount (default: 0)
    --hqStartTime HQSTARTTIME
                          ZMW Filter: exclude ZMW if start time of HQ
                          regiongreater than this value (seconds) (default: inf)
    --minReadScore MINREADSCORE
                          ZMW Filter: exclude ZMW if readScore is less thanthis
                          value (default: 0)
    --minAvgBarcodeScore MINAVGBARCODESCORE
                          ZMW Filter: exclude ZMW if average barcode score is
                          less than this value (default: 0.0)
    --minNumBarcodes MINNUMBARCODES
                          ZMW Filter: exclude ZMW if number of barcodes observed
                          is less than this value (default: 1)
    --minScoreRatio MINSCORERATIO
                          ZMW Filter: exclude ZMWs whose best score divided by
                          the 2nd best score is less than this ratio (default:
                          1.0)
                          

The ``emitFastqs`` command takes as input both an input.fofn for the
bas.h5 files as well as a barcode.fofn from a call to labelZmws. The
optional parameter ``outDir`` dictates where the files will be
written. For each detected barcode, a fast[a|q] file will be emitted
with all of the reads for that barcode. The ``trim`` parameter
dictates how much of the read should be trimmed off. The default
parameter for ``trim`` is the length of the barcode (which is stored
in the barcode hdf5 files). At the moment, all barcodes in the barcode
FASTA file must be the same length, therefore only a constant trim
value is supported. In practice, one can aggressively trim in order to
ensure that extra bases aren't left on the ends of reads. Finally, the
``subreads`` parameter dictates whether subreads or CCS reads should
be returned with the default being the appropriate reads according to
the input file type, either CCS or subreads. This parameter is only
inspected if the input.fofn contains both CCS and subread data, if the
input.fofn contains only subread or CCS data then that is returned
irrespective of the state of the the ``subreads`` parameter and a
warning is issued.

consensus
---------
  usage: pbbarcode.py consensus [-h] [--subsample SUBSAMPLE] [--nZmws NZMWS]
                                [--outDir OUTDIR] [--keepTmpDir]
                                [--ccsFofn CCSFOFN] [--nProcs NPROCS]
                                [--noQuiver]
                                [--minMaxInsertLength MINMAXINSERTLENGTH]
                                [--hqStartTime HQSTARTTIME]
                                [--minReadScore MINREADSCORE]
                                [--minAvgBarcodeScore MINAVGBARCODESCORE]
                                [--minNumBarcodes MINNUMBARCODES]
                                [--minScoreRatio MINSCORERATIO]
                                [--barcode BARCODE [BARCODE ...]]
                                input.fofn barcode.fofn

  Compute consensus sequences for each barcode.

  positional arguments:
    input.fofn            input bas.h5 fofn file
    barcode.fofn          input bc.h5 fofn file

  optional arguments:
    -h, --help            show this help message and exit
    --subsample SUBSAMPLE
                          Subsample ZMWs (default: 1)
    --nZmws NZMWS         Take n ZMWs (default: -1)
    --outDir OUTDIR       Use this directory to output results (default: .)
    --keepTmpDir
    --ccsFofn CCSFOFN     Obtain CCS data from ccsFofn instead of input.fofn
                          (default: )
    --nProcs NPROCS       Use nProcs to execute. (default: 16)
    --noQuiver
    --minMaxInsertLength MINMAXINSERTLENGTH
                          ZMW Filter: exclude ZMW if the longest subreadis less
                          than this amount (default: 0)
    --hqStartTime HQSTARTTIME
                          ZMW Filter: exclude ZMW if start time of HQ
                          regiongreater than this value (seconds) (default: inf)
    --minReadScore MINREADSCORE
                          ZMW Filter: exclude ZMW if readScore is less thanthis
                          value (default: 0)
    --minAvgBarcodeScore MINAVGBARCODESCORE
                          ZMW Filter: exclude ZMW if average barcode score is
                          less than this value (default: 0.0)
    --minNumBarcodes MINNUMBARCODES
                          ZMW Filter: exclude ZMW if number of barcodes observed
                          is less than this value (default: 1)
    --minScoreRatio MINSCORERATIO
                          ZMW Filter: exclude ZMWs whose best score divided by
                          the 2nd best score is less than this ratio (default:
                          1.0)
    --barcode BARCODE [BARCODE ...]
                          Use this to extract consensus for just one barcode.
                          (default: None)

The ``emitFastqs`` command takes as input both an input.fofn for the
bas.h5 files as well as a barcode.fofn from a call to labelZmws. The
results are a FASTA file with an entry for each barcode containing the
consensus amplicon sequence. This mode utilizes ``Quiver`` and
``pbdagcon`` to compute consensus.   

In cases where the amplicon is fewer than 2.5k bases, using CCS data
is quite helpful. The ``--ccsFofn`` allows one to pass directly the
ccs files. In many cases, both the CCS and raw basecalls are in the
same file so you can check by passing the same parameter to input.fofn
as to ccsFofn. 

Dependencies
````````````

The pbbarcode package depends on a standard pbcore installation
(https://github.com/PacificBiosciences/pbcore). If one wishes to use
the ``consensus`` tool, ``pbdagcon`` needs to be installed
(https://github.com/PacificBiosciences/pbdagcon).


Barcode HDF5 File
`````````````````

The barcode hdf5 file, ``bc.h5``, represents a simple data store for
barcode calls and their scores for each ZMW. Generally, a user need
not interact with barcode hdf5 files, but can use the results stored in
either the resulting cmp.h5 file or fast[a|q] files. The barcode hdf5
file contains the following structure:

/BarcodeCalls/best - (nZMWs, 6)[32-bit integer] dataset with the
following columns: 

    ``holeNumber,nAdapters,barcodeIdx1,barcodeScore1,barcodeIdx2,barcodeScore2``

Additionally, the ``best`` dataset has the following attributes:

+-----------+------------------------------------------------------------------------+
|movieName  |m120408_042614_richard_c100309392550000001523011508061222_s1_p0         |
|           |                                                                        |
+-----------+------------------------------------------------------------------------+
|columnNames|holeNumber,nAdapters,barcodeIdx1,barcodeScore1,barcodeIdx2,             |
|           |barcodeScore2                                                           |
+-----------+------------------------------------------------------------------------+
|scoreMode  |[symmetric|paired]                                                      |
|           |                                                                        |
+-----------+------------------------------------------------------------------------+
|barcodes   |'bc_1', 'bc_2', ...., 'bc_N'                                            |
|           |                                                                        |
+-----------+------------------------------------------------------------------------+

The two barcodeIdx1 and barcodeIdx2 columns are indices into
``barcodes`` attribute. The ``scoreMode`` is scoring mode used to
align the barcodes. The ``barcodes`` attribute correspond to the
barcode.fasta sequence names. 

Additionally, in some circumstances, it is useful to retain the entire
history of the scoring, i.e., each barcode scored to each adapter
across all ZMWs. In oder to retain this information, one must call:

    ``pbbarcode.py labelZmws --saveExtendedInfo ...``

In this mode, the resultant HDF5 file will have an additional dataset
under the BarcodeCalls group, named: ``all``. This dataset has the
following format:

/BarcodeCalls/all - (nbarcodes * nadapters[zmw_i], 4) \forall i in 1 ... nZMWs 

    ```holeNumber, adapterIdx, barcodeIdx, score```

The ``adapterIdx`` is the index of the adapter along the molecule,
i.e., adapterIdx 1 is the first adapter scored.

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
