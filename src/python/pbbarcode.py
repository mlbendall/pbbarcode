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

import h5py as h5
import numpy as n

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.BasH5Reader import *
from pbcore.io.CmpH5Reader import *
from pbtools.pbbarcode.BarcodeLabeler import BarcodeScorer
from pbcore.io.FastaIO import *

__version__ = ".04"

## Paths to the Barcode Datasets.
BC_DS_PATH     = "BarcodeCalls/best"
BC_INFO_NAME   = "BarcodeInfo/Name"
BC_INFO_ID     = "BarcodeInfo/ID"
BC_ALN_INFO_DS = "AlnInfo/Barcode"

NULL_BARCODE   = "<>"

def makeBCLabel(s1, s2):
    return NULL_BARCODE.join((s1, s2))
    
class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))
        subparsers = self.getSubParsers()
        
        desc = ['Creates a barcode.h5 file from base h5 files.']
        parser_m = subparsers.add_parser('labelZMWs', description = "\n".join(desc), parents = [self.parser])
        parser_m.add_argument('--outDir', help = 'Where to write the newly created barcode.h5 files.',
                              default = os.getcwd())
        parser_m.add_argument('--outFofn', help = 'Write to outFofn',
                              default = 'barcode.fofn')
        parser_m.add_argument('--adapterSidePad', help = 'Pad with adapterSidePad bases',
                              default = 0, type = int)
        parser_m.add_argument('--insertSidePad', help = 'Pad with insertSidePad bases',
                              default = 4, type = int)
        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta', 
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.fofn',
                              help = 'Input base fofn')
        
        desc = ['adds information about barcode alignments to a cmp.h5 file',
                'from a previous call to "labelZMWs".']
        parser_s = subparsers.add_parser('labelAlignments', description = "\n".join(desc),
                                         parents = [self.parser], help = "Label reads from a barcode or region h5 file")
        parser_s.add_argument('--asymmetric', help = 'whether the same barcode is on both sides of the molecule.',
                              action="store_true", default = False)
        parser_s.add_argument('inputFofn', metavar = 'barcode.fofn',
                              help = 'input barcode fofn file')
        parser_s.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                              help = 'cmp.h5 file to add barcode labels')
        
    def getVersion(self):
        return __version__

    def makeBarcodeH5FromBasH5(self, basH5):
        oFile = '/'.join((self.args.outDir, basH5.movieName + '.bc.h5'))
        logging.info("Labeling regions, preparing to write: %s" % oFile) 
        labeler = BarcodeScorer(basH5, FastaReader(self.args.barcodeFile),
                                self.args.adapterSidePad, self.args.insertSidePad)
        allScores = labeler.scoreZMWs()
        logging.debug("Scored ZMWs")
        bestScores = labeler.chooseBestBarcodes(allScores)

        ## write the new file. 
        outH5  = h5.File(oFile, 'a')
        outDta = n.vstack(bestScores)
        if BC_DS_PATH in outH5:
            del outH5[BC_DS_PATH]
        bestDS = outH5.create_dataset(BC_DS_PATH, data = outDta, dtype = "int32")
        bestDS.attrs['movieName'] = basH5.movieName
        bestDS.attrs['barcodes'] = n.array([x.name for x in FastaReader(self.args.barcodeFile)], 
                                           dtype = h5.new_vlen(str))
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
        def movieName(movieFile):
            return os.path.basename(movieFile).replace('.bc.h5', '')
        
        def makeMovieMap(movieFile):
             bcFile = h5.File(movieFile, 'r')
             calls = bcFile[BC_DS_PATH][:]
             bcFile.close()
             return dict(zip(calls[:,0], calls[:,1:calls.shape[1]]))
            
        def getBarcodeLabels(movieFile):
            bcFile = h5.File(movieFile, 'r')
            labels = bcFile[BC_DS_PATH].attrs['barcodes'][:]
            bcFile.close()
            return labels
            
        logging.info("Labeling alignments using %s labeling." % 
                     'asymmetric' if self.args.asymmetric else 'symmetric')
    
        barcodeFofn = open(self.args.inputFofn).read().splitlines()
        movieMap = {movieName(movie) : makeMovieMap(movie) for movie in barcodeFofn}
        cmpH5 = CmpH5Reader(self.args.cmpH5)
        labels = getBarcodeLabels(barcodeFofn[0])
        bcDS = n.zeros((len(cmpH5), 3), dtype = "int32")

        if not self.args.asymmetric:
            bcLabels = [ makeBCLabel(a, b) for a,b in zip(labels, labels)]
            bcLabels.append(NULL_BARCODE)
            for (i, aln) in enumerate(cmpH5):
                try:
                    mp = movieMap[aln.movieInfo.Name][aln.HoleNumber]
                    bcDS[i,0] = mp[1] # index
                    bcDS[i,1] = mp[2] # score
                    bcDS[i,2] = mp[0] # count
                except:
                    bcDS[i,0] = len(bcLabels) - 1 # index of NULL_BARCODE. 
                    bcDS[i,1] = 0
                    bcDS[i,2] = 0
        else:
            bcLabels = [makeBCLabel(labels[i], labels[j]) for i in xrange(0, len(labels)) 
                        for j in xrange(i+1, len(labels))]
            bcLabels.append(NULL_BARCODE)
            for (i, aln) in enumerate(cmpH5):
                try:
                    mp = movieMap[aln.movieInfo.Name][aln.HoleNumber]
                    # this is a little trickier. 
                    (l, h) = (mp[1], mp[3]) if mp[1] < mp[3] else (mp[3], mp[1])
                    idx = ((len(labels)-1)*l - l*(l-1)/2) + h-l-1
                    bcDS[i,0] = idx         # index
                    bcDS[i,1] = mp[2]+mp[4] # score
                    bcDS[i,2] = mp[0]       # count
                    if not makeBCLabel(labels[l], labels[h]) == bcLabels[idx]:
                        raise Exception("Software Error - not equal: %s to %s" % (makeBCLabel(labels[l], labels[h]), 
                                                                                  bcLabels[idx]))
                except:
                    bcDS[i,0] = len(bcLabels) - 1 # index of NULL_BARCODE.
                    bcDS[i,1] = 0
                    bcDS[i,2] = 0
            
        # Now, we write the datastructures to the file.
        try:
            cmpH5.__del__()
        except Exception, e:
            log.info(e) # this shouldn't happen.

        H5 = h5.File(self.args.cmpH5, 'r+')
        if BC_INFO_ID in H5:
            del H5[BC_INFO_ID]
        if BC_INFO_NAME in H5:
            del H5[BC_INFO_NAME]
        H5.create_dataset(BC_INFO_ID, data = n.array(range(0, len(bcLabels))), dtype = 'int32')
        H5.create_dataset(BC_INFO_NAME, data = bcLabels, dtype = h5.new_vlen(str))
        
        if BC_ALN_INFO_DS in H5:
            del H5[BC_ALN_INFO_DS]
        bcDS = H5.create_dataset(BC_ALN_INFO_DS, data = bcDS, dtype = 'int32')
        bcDS.attrs['ColumnNames'] = n.array(['index', 'score', 'count'])
        bcDS.attrs['BarcodeMode'] = "asymmetric" if self.args.asymmetric else "symmetric"
        H5.close()

    def run(self):
        logging.debug("Arguments" + str(self.args))
        if self.args.subName == 'labelZMWs':
            self.makeBarcodeFofnFromBasFofn()
        elif self.args.subName == 'labelAlignments':
            self.labelAlignments()

if __name__ == '__main__':    
    sys.exit(Pbbarcode().start())
    


        

