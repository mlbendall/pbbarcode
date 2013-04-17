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

#!/usr/bin/env python
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
from pbtools.pbbarcode.BarcodeH5Reader import BarcodeH5Reader, BC_DS_PATH, \
    BarcodeIdxException
from pbcore.io import FastaReader, FastqWriter, FastqRecord

__version__ = ".08"

## Paths to the Barcode Datasets.
BC_INFO_NAME   = "BarcodeInfo/Name"
BC_INFO_ID     = "BarcodeInfo/ID"
BC_ALN_INFO_DS = "AlnInfo/Barcode"
SCORE_MODES    = ['symmetric', 'paired', 'asymmetric']

class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))
        subparsers = self.subParsers
        
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
        parser_m.add_argument('--nZMWs', type = int, default = -1, help = 'Use nZMWs for testing')
        parser_m.add_argument('--maxAdapters', type = int, default = 20, 
                              help = 'If there are more than maxAdapters, ignore')
        parser_m.add_argument('--scoreFirst', action = 'store_true', default = False,
                              help = 'Whether to try to score the leftmost barcode in a trace.')

        parser_m.add_argument('--new', action = 'store_true', default = False)

        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta', 
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.fofn',
                              help = 'Input base fofn')
        
        desc = ['adds information about barcode alignments to a cmp.h5 file',
                'from a previous call to "labelZMWs".']
        parser_s = subparsers.add_parser('labelAlignments', description = "\n".join(desc),
                                         help = "Label reads from a barcode or region h5 file")
        parser_s.add_argument('inputFofn', metavar = 'barcode.fofn',
                              help = 'input barcode fofn file')
        parser_s.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                              help = 'cmp.h5 file to add barcode labels')
        
        desc = ['Takes a bas.h5 fofn and a barcode.h5 fofn and produces',
                'a fastq file for each barcode found.']
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
                              action = "store_true",
                              default = False)

        desc = ['Takes a bas.h5 fofn and a barcode.h5 fofn and produces',
                'a fasta file for each barcode with a consensus sequence.']
        parser_s = subparsers.add_parser('consensus', description = "\n".join(desc),
                                         help = "compute consensus")
        parser_s.add_argument('inputFofn', metavar = 'input.fofn',
                              help = 'input bas.h5 fofn file')
        parser_s.add_argument('barcodeFofn', metavar = 'barcode.fofn',
                              help = 'input barcode.h5 fofn file')
        parser_s.add_argument('--trim', help = 'trim off barcodes and any excess constant sequence',
                              default = 20, type = int)
        parser_s.add_argument('--subreads', help = ('whether to produce fastq files for the subreads;' +
                                                    'the default is to use the CCS reads.'),
                              action = "store_true",
                              default = False)

    def getVersion(self):
        return __version__

    def makeBarcodeH5FromBasH5(self, basH5):
        # check to make sure we can write the output before we do a
        # bunch of work.
        oFile = '/'.join((self.args.outDir, re.sub(r'\.ba[x|s]\.h5', '.bc.h5', 
                                                   os.path.basename(basH5.filename))))
        logging.debug("Setting up output file: %s" % oFile)
        outH5 = h5.File(oFile, 'a')

        logging.info("Labeling: %s" % basH5.filename)
        labeler = BarcodeScorer(basH5, FastaReader(self.args.barcodeFile),
                                self.args.adapterSidePad, self.args.insertSidePad,
                                scoreMode = self.args.scoreMode, maxHits = self.args.maxAdapters,
                                scoreFirst = self.args.scoreFirst, startTimeCutoff = 2, 
                                scorePairedNew = self.args.new)
        
        if self.args.nZMWs < 0:
            zmws = basH5.sequencingZmws
        else:
            zmws = basH5.sequencingZmws[0:self.args.nZMWs]

        allScores = labeler.scoreZMWs(zmws = zmws)
        logging.debug("Scored ZMWs")

        bestScores = labeler.chooseBestBarcodes(allScores)
        logging.debug("Chose best barcodes")

        ## write the new file. 
        logging.info("Writing to: %s" % oFile)
        outDta = n.vstack(bestScores)
        if BC_DS_PATH in outH5:
            del outH5[BC_DS_PATH]
        bestDS = outH5.create_dataset(BC_DS_PATH, data = outDta, dtype = "int32")
        bestDS.attrs['movieName'] = basH5.movieName
        bestDS.attrs['barcodes'] = n.array([x.name for x in FastaReader(self.args.barcodeFile)], 
                                           dtype = h5.new_vlen(str))
        bestDS.attrs['columnNames'] = n.array(['holeNumber', 'nAdapters', 'barcodeIdx1', 
                                               'barcodeScore1', 'barcodeIdx2', 'barcodeScore2'], 
                                              dtype = h5.new_vlen(str))
        bestDS.attrs['scoreMode'] = self.args.scoreMode


        outH5.close()
        return oFile

    def makeBarcodeFofnFromBasFofn(self):
        inputFofn = self.args.inputFile
        inFiles = open(inputFofn).read().splitlines()
        newFiles = [self.makeBarcodeH5FromBasH5(BasH5Reader(basH5)) for basH5 in inFiles]
        oFile = open(self.args.outFofn, 'w')
        for nF in newFiles:
            oFile.write(nF + "\n")
        oFile.close()

    def labelAlignments(self):
        logging.info("Labeling alignments using: %s" % self.args.inputFofn)
        def movieName(movieFile):
            r = os.path.basename(movieFile).replace('.bc.h5', '')
            # ibid, bug 23071
            return re.sub(r'\.[1-9]$', '', r)
            
        barcodeFofn = open(self.args.inputFofn).read().splitlines()
        movieMap = { movieName(movie) : BarcodeH5Reader(movie) for movie in barcodeFofn }

        with CmpH5Reader(self.args.cmpH5) as cmpH5:
            bcDS = n.zeros((len(cmpH5), 3), dtype = "int32")
            for (i, aln) in enumerate(cmpH5):
                bcFile = movieMap[aln.movieInfo.Name]
                zmwCall = bcFile.getBarcodeTupleForZMW(aln.HoleNumber)
                bcDS[i,:] = n.array(zmwCall)

        H5 = h5.File(self.args.cmpH5, 'r+')
        if BC_INFO_ID in H5:
            del H5[BC_INFO_ID]
        if BC_INFO_NAME in H5:
            del H5[BC_INFO_NAME]
     
        # we use the first one to get the labels, if somehow they
        # don't have all of the same stuff that will be an issue.
        bcReader = movieMap[movieMap.keys()[0]]
        bcLabels = bcReader.bcLabels
        H5.create_dataset(BC_INFO_ID, data = n.array(range(0, len(bcLabels))), dtype = 'int32')
        H5.create_dataset(BC_INFO_NAME, data = bcLabels, dtype = h5.new_vlen(str))
        if BC_ALN_INFO_DS in H5:
            del H5[BC_ALN_INFO_DS]
        bcDS = H5.create_dataset(BC_ALN_INFO_DS, data = bcDS, dtype = 'int32')
        bcDS.attrs['ColumnNames'] = n.array(['index', 'score', 'count'])
        bcDS.attrs['BarcodeMode'] = bcReader.scoreMode
        H5.close()
        

    def getFastqs(self):
        # step through the bas.h5 and barcode.h5 files and emit 
        # reads for each of these. 
        inputFofn = n.array(open(self.args.inputFofn).read().splitlines())
        barcodeFofn = n.array(open(self.args.barcodeFofn).read().splitlines())

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
            bcH5  = BarcodeH5Reader(barcodeFile)
            
            for label in bcH5.bcLabels:
                try:
                    zmws = bcH5.getZMWsForBarcode(label)
                    for row in range(0, zmws.shape[0]):
                        zmw = basH5[zmws[row, 0]]
                        if zmw:
                            if self.args.subreads:
                                reads = zmw.subreads()
                            else:
                                reads = [zmw.ccsRead]
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

    def emitFastqs(self):
        outFiles = self.getFastqs()
        outDir   = self.args.outDir
        if not os.path.exists(self.args.outDir):
            os.makedirs(self.args.outDir)
        
        for k in outFiles.keys():
            with FastqWriter("%s/%s.fastq" % (self.args.outDir, k)) as w:
                for e in outFiles[k]:
                    ## XXX : Empty Strings!! 
                    w.writeRecord(trimFastqRecord(e, self.args.trim))

    def run(self):
        logging.debug("Arguments" + str(self.args))
        if self.args.subCommand == 'labelZMWs':
            self.makeBarcodeFofnFromBasFofn()
        elif self.args.subCommand == 'labelAlignments':
            self.labelAlignments()
        elif self.args.subCommand == 'emitFastqs':
            self.emitFastqs()
        elif self.args.subCommand == 'consensus':
            self.consensus()
        else:
            sys.exit(1)

def trimFastqRecord(fastqRecord, trim):
    return FastqRecord(fastqRecord.name,
                       fastqRecord.sequence[trim:-trim],
                       fastqRecord.quality[trim:-trim])

            
if __name__ == '__main__':    
    sys.exit(Pbbarcode().start())
    


        

