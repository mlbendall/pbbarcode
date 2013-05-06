#!/usr/bin/env python
#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os
import sys
import argparse
import logging
import tempfile
import shutil
import pkg_resources
import re

import h5py as h5
import numpy as n

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.BasH5Reader import *
from pbcore.io.CmpH5Reader import *
from pbtools.pbbarcode.BarcodeLabeler import BarcodeScorer
from pbtools.pbbarcode.BarcodeH5Reader import BC_DS_PATH, BarcodeIdxException,\
    BC_DS_ALL_PATH, create, BarcodeH5Reader
from pbcore.io import FastaReader, FastqWriter, FastqRecord, \
    FastaWriter, FastaRecord

__version__ = ".1"

## Paths to the Barcode Datasets.
BC_INFO_NAME   = "BarcodeInfo/Name"
BC_INFO_ID     = "BarcodeInfo/ID"
BC_ALN_INFO_DS = "AlnInfo/Barcode"
SCORE_MODES    = ['symmetric', 'paired', 'asymmetric']

def makeBarcodeH5FromBasH5(basH5):
    # check to make sure we can write the output before we do a
    # bunch of work.
    oFile = '/'.join((runner.args.outDir, re.sub(r'\.ba[x|s]\.h5', '.bc.h5', 
                                               os.path.basename(basH5.filename))))
    logging.debug("Setting up output file: %s" % oFile)
    outH5 = h5.File(oFile, 'a')

    logging.info("Labeling: %s" % basH5.filename)
    labeler = BarcodeScorer(basH5, FastaReader(runner.args.barcodeFile),
                            runner.args.adapterSidePad, runner.args.insertSidePad,
                            scoreMode = runner.args.scoreMode, 
                            maxHits = runner.args.maxAdapters,
                            scoreFirst = runner.args.scoreFirst, 
                            startTimeCutoff = 1) 
    if runner.args.nZMWs < 0:
        zmws = basH5.sequencingZmws
    else:
        zmws = basH5.sequencingZmws[0:runner.args.nZMWs]

    allScores = labeler.scoreZMWs(zmws = zmws)
    logging.debug("Scored ZMWs")

    bestScores = labeler.chooseBestBarcodes(allScores)
    logging.debug("Chose best barcodes")

    ## write the new file. 
    logging.info("Writing to: %s" % oFile)

    ## best data set chosen according to the scoring mode.
    outDta = n.vstack(bestScores)
    if BC_DS_PATH in outH5:
        del outH5[BC_DS_PATH]
    bestDS = outH5.create_dataset(BC_DS_PATH, data = outDta, dtype = "int32")
    bestDS.attrs['movieName'] = basH5.movieName
    bestDS.attrs['barcodes'] = n.array([x.name for x in FastaReader(runner.args.barcodeFile)], 
                                       dtype = h5.new_vlen(str))
    bestDS.attrs['columnNames'] = n.array(['holeNumber', 'nAdapters', 'barcodeIdx1', 
                                           'barcodeScore1', 'barcodeIdx2', 'barcodeScore2'], 
                                          dtype = h5.new_vlen(str))
    bestDS.attrs['scoreMode'] = runner.args.scoreMode

    ## all dataset.
    def makeRecord(r):
        #  (zmw, adapter, barcode, score)
        def makeArray(l, v):
            a = n.zeros(l, dtype = type(v))
            a[:] = v
            return a
        z = r[0]
        c = r[3]
        c = [ i.round() for i in c ]
        nbarcodes = len(c[0])
        nadapters = len(c)
        zmws = makeArray(nbarcodes*nadapters, z)
        adapters = n.concatenate([ makeArray(nbarcodes, i) for i in xrange(1, nadapters + 1)])
        idxs = n.concatenate([ range(0, nbarcodes) for i in xrange(0, nadapters)])
        scores = n.concatenate(c)
        return n.transpose(n.vstack((zmws, adapters, idxs, scores)))

    records = [ makeRecord(r) for r in allScores if r[1] ] # if there are any scores to look at.
    records = n.vstack(records)
    
    if BC_DS_ALL_PATH in outH5:
        del outH5[BC_DS_ALL_PATH]
    allDS = outH5.create_dataset(BC_DS_ALL_PATH, data = records, dtype = 'int32')
    allDS.attrs['movieName'] = basH5.movieName
    allDS.attrs['barcodes'] = n.array([x.name for x in FastaReader(runner.args.barcodeFile)], 
                                      dtype = h5.new_vlen(str))
    allDS.attrs['columnNames'] = n.array(['holeNumber', 'adapter', 'barcodeIdx', 'score'], 
                                         dtype = h5.new_vlen(str))

    outH5.close()
    return oFile

def mpWrapper(f):
    return makeBarcodeH5FromBasH5(BasH5Reader(f))

def makeBarcodeFofnFromBasFofn():
    inputFofn = runner.args.inputFile
    inFiles = open(inputFofn).read().splitlines()
    if runner.args.nProcs <= 1:
        newFiles = map(mpWrapper, inFiles)
    else:
  #      pool = Pool(runner.args.nProcs)

#        logging.debug("Using %d processes." % runner.args.nProcs)
        newFiles = map(mpWrapper, inFiles)

    oFile = open(runner.args.outFofn, 'w')
    for nF in newFiles:
        oFile.write(nF + "\n")
    oFile.close()

def labelAlignments():
    logging.info("Labeling alignments using: %s" % runner.args.inputFofn)

    def movieBasePath(bcFile):
        r = os.path.basename(bcFile).replace('.bc.h5', '')
        return '/'.join((os.path.dirname(bcFile), re.sub(r'\.[1-9]$', '', r)))

    barcodeFofn = open(runner.args.inputFofn).read().splitlines()
    movieNames  = n.unique(map(movieBasePath, barcodeFofn))
    movieMap    = { os.path.basename(movie):create(movie) for movie in movieNames }

    with CmpH5Reader(runner.args.cmpH5) as cmpH5:
        bcDS = n.zeros((len(cmpH5), 3), dtype = "int32")
        for (i, aln) in enumerate(cmpH5):
            bcFile = movieMap[aln.movieInfo.Name]
            zmwCall = bcFile.getBarcodeTupleForZMW(aln.HoleNumber)
            bcDS[i,:] = n.array(zmwCall)

    H5 = h5.File(runner.args.cmpH5, 'r+')
    if BC_INFO_ID in H5:
        del H5[BC_INFO_ID]
    if BC_INFO_NAME in H5:
        del H5[BC_INFO_NAME]

    # we use the first one to get the labels, if somehow they
    # don't have all of the same stuff that will be an issue.
    bcReader = movieMap[movieMap.keys()[0]]
    bcLabels = bcReader.getBarcodeLabels()
    H5.create_dataset(BC_INFO_ID, data = n.array(range(0, len(bcLabels))), dtype = 'int32')
    H5.create_dataset(BC_INFO_NAME, data = bcLabels, dtype = h5.new_vlen(str))
    if BC_ALN_INFO_DS in H5:
        del H5[BC_ALN_INFO_DS]
    bcDS = H5.create_dataset(BC_ALN_INFO_DS, data = bcDS, dtype = 'int32')
    bcDS.attrs['ColumnNames'] = n.array(['index', 'score', 'count'])
    bcDS.attrs['BarcodeMode'] = bcReader.getScoreMode()
    H5.close()


def getFastqs():
    # step through the bas.h5 and barcode.h5 files and emit 
    # reads for each of these. 
    inputFofn = n.array(open(runner.args.inputFofn).read().splitlines())
    barcodeFofn = n.array(open(runner.args.barcodeFofn).read().splitlines())

    def movieNameFromFile(fn): 
        return re.sub(r'\.bc\.h5|\.bax\.h5|\.bas\.h5', '', os.path.basename(fn))

    # sort them.
    inputFofn = list(inputFofn[n.array(n.argsort([movieNameFromFile(a) for \
                                                      a in inputFofn]))])
    barcodeFofn = list(barcodeFofn[n.array(n.argsort([movieNameFromFile(a) for \
                                                          a in barcodeFofn]))])

    if len(inputFofn) != len(barcodeFofn) or \
            any([ movieNameFromFile(a) != movieNameFromFile(b) \
                      for a,b in zip(inputFofn, barcodeFofn)]):
        raise Exception("input.fofn and barcode.fofn do not match.")

    outFiles = {} 
    for basFile, barcodeFile in zip(inputFofn, barcodeFofn):
        logging.info("Processing: %s %s" % (basFile, barcodeFile))
        basH5 = BasH5Reader(basFile)
        # XXX: here, I'll use this interface because it is easier
        bcH5  = BarcodeH5Reader(barcodeFile)

        for label in bcH5.bcLabels:
            try:
                zmws = bcH5.getZMWsForBarcode(label)
                for row in range(0, zmws.shape[0]):
                    zmw = basH5[zmws[row, 0]]
                    if zmw:
                        if runner.args.subreads:
                            reads = zmw.subreads()
                        else:
                            reads = [zmw.ccsRead()]
                        if any(reads):
                            if not label in outFiles.keys():
                                outFiles[label] = []
                            for read in reads:
                                outFiles[label].append(FastqRecord(read.readName, 
                                                                   read.basecalls(),
                                                                   read.QualityValue()))
            except BarcodeIdxException, e:
                continue

    return outFiles

def emitFastqs():
    outFiles = getFastqs()
    outDir   = runner.args.outDir
    fasta    = runner.args.fasta

    if not os.path.exists(runner.args.outDir):
        os.makedirs(runner.args.outDir)

    if fasta:
        writer = FastaWriter
        def record(n, s, qv):
            FastaRecord(n, s)
    else:
        writer = FastqWriter
        record = FastqRecord

    for k in outFiles.keys():
        with writer("%s/%s.fastq" % (runner.args.outDir, k)) as w:
            for e in outFiles[k]:
                r = record(e.name,
                           e.sequence[runner.args.trim:-runner.args.trim],
                           e.quality[runner.args.trim:-runner.args.trim])
                if r:
                    w.writeRecord(r)


class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))
        subparsers = self.getSubParsers()
        
        desc = ['Creates a barcode.h5 file from base h5 files.']
        parser_m = subparsers.add_parser('labelZMWs', description = "\n".join(desc), 
                                         help = 'Label zmws with barcode annotation')
        parser_m.add_argument('--outDir', help = 'Where to write the newly created barcode.h5 files.',
                              default = os.getcwd())
        parser_m.add_argument('--outFofn', help = 'Write to outFofn',
                              default = 'barcode.fofn')
        parser_m.add_argument('--adapterSidePad', help = 'Pad with adapterSidePad bases',
                              default = 4, type = int)
        parser_m.add_argument('--insertSidePad', help = 'Pad with insertSidePad bases',
                              default = 4, type = int)
        parser_m.add_argument('--scoreMode', help = 'The mode in which the barcodes should be scored.',
                              choices = SCORE_MODES, default = 'symmetric', 
                              type = str)
        parser_m.add_argument('--maxAdapters', type = int, default = 20, 
                              help = 'Only score the first maxAdapters')
        parser_m.add_argument('--scoreFirst', action = 'store_true', default = False,
                              help = 'Whether to try to score the leftmost barcode in a trace.')
        parser_m.add_argument('--nZMWs', type = int, default = -1, help = 'Use nZMWs for testing')
        parser_m.add_argument('--nProcs', type = int, default = 8, help = 'How many processes to use')
        
        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta', 
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.fofn',
                              help = 'Input base fofn')
        
        desc = ['Adds information about barcode alignments to a cmp.h5 file',
                'from a previous call to "labelZMWs".']
        parser_s = subparsers.add_parser('labelAlignments', description = "\n".join(desc),
                                         help = "Label reads from a barcode or region h5 file")
        parser_s.add_argument('inputFofn', metavar = 'barcode.fofn',
                              help = 'input barcode fofn file')
        parser_s.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                              help = 'cmp.h5 file to add barcode labels')
        
        desc = ['Takes a bas.h5 fofn and a barcode.h5 fofn and produces',
                'a fast[a|q] file for each barcode.']
        parser_s = subparsers.add_parser('emitFastqs', description = "\n".join(desc),
                                         help = "Write fastq files")
        parser_s.add_argument('inputFofn', metavar = 'input.fofn',
                              help = 'input bas.h5 fofn file')
        parser_s.add_argument('barcodeFofn', metavar = 'barcode.fofn',
                              help = 'input barcode.h5 fofn file')
        parser_s.add_argument('--outDir', metavar = 'output.dir',
                              help = 'output directory to write fastq files',
                              default = os.getcwd())
        parser_s.add_argument('--trim', help = 'trim off barcodes and any excess constant sequence',
                              default = 20, type = int)
        parser_s.add_argument('--subreads', help = ('whether to produce fastq files for the subreads;' +
                                                    'the default is to use the CCS reads.'),
                              action = 'store_true',
                              default = False)
        parser_s.add_argument('--fasta', help = ('whether the files produced should be FASTA files as' +
                                                 'opposed to FASTQ'),
                              action = 'store_true',
                              default = False)

    def getVersion(self):
        return __version__

    def run(self):
        logging.debug("Arguments" + str(self.args))
        
        if self.args.subName == 'labelZMWs':
            makeBarcodeFofnFromBasFofn()
        elif self.args.subName == 'labelAlignments':
            labelAlignments()
        elif self.args.subName == 'emitFastqs':
            emitFastqs()
        else:
            sys.exit(1)

           
if __name__ == '__main__':    
    # runner becomes a global to facilitate 'multiprocessing'
    runner = Pbbarcode()
    sys.exit(runner.start())
    


        

