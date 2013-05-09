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
import logging

import h5py as h5
import numpy as n

from pbcore.io.BasH5Reader import *
from pbcore.io.CmpH5Reader import *

NULL_BARCODE   = "--"
BC_DS_PATH     = "BarcodeCalls/best"
BC_DS_ALL_PATH = "BarcodeCalls/all"

def makeBCLabel(s1, s2):
    return NULL_BARCODE.join((s1, s2))

class BarcodeIdxException(Exception):
    pass

def create(movieName):
    ## XXX : Detecting multi-part; pretty ugly. 
    if os.path.exists(movieName + '.bc.h5'):
        logging.debug("Instantiating BarcodeH5Reader")
        return BarcodeH5Reader(movieName + '.bc.h5')
    else:
        logging.debug("Instantiating MPBarcodeReader")
        parts = map(lambda z : '.'.join((movieName, str(z), 'bc.h5')), [1,2,3])
        logging.debug("Trying to load parts:" + '\n'.join(parts))
        parts = filter(lambda p : os.path.exists(p), parts)

        if parts:
            return MPBarcodeH5Reader(parts)
        else:
            raise Exception("Unable to instantiate a BarcodeH5Reader")

class MPBarcodeH5Reader(object):
    def __init__(self, fnames):
        self.parts = map(BarcodeH5Reader, fnames)
        def rng(x):
            return (n.min(x), n.max(x))
        # these aren't the ranges of ZMWs, but the ranges for the
        # scored ZMWs.
        self.bins = map(lambda z : rng(z.bestDS[:,0]), self.parts)

    def choosePart(self, holeNumber):
        for i,b in enumerate(self.bins):
            if holeNumber >= b[0] and holeNumber <= b[1]:
                return self.parts[i]
        # Return None meaning the zmw is ouf of the range of
        # the scored ZMWs for all parts.
        return None

    def getBarcodeTupleForZMW(self, holeNumber):
        part = self.choosePart(holeNumber)
        if part:
            return part.getBarcodeTupleForZMW(holeNumber)
        else:
            # I have not scored this alignment, therefore we return the
            # NULL tuple
            return self.parts[0].nullBarcodeTuple()

    def getExtendedBarcodeInfoForZMW(self, holeNumber):
        part = self.choosePart(holeNumber)
        if part:
            return part.getExtendedBarcodeInfoForZMW(holeNumber)
        else:
            return None

    def getZMWsForBarcode(self, barcodeName):
        raise NotImplementedError("Directly use BarcodeH5Reader for this task")
    
    def getBarcodeLabels(self):
        return self.parts[0].bcLabels

    def getScoreMode(self):
        return self.parts[0].scoreMode
        
class BarcodeH5Reader(object):
    def __init__(self, fname):
        self.h5File    = h5.File(fname, 'r')
        self.bestDS    = self.h5File[BC_DS_PATH]
        self.scoreMode = self.bestDS.attrs['scoreMode']
        self.barcodes  = self.bestDS.attrs['barcodes']
        ## note: I hack off the ZMW.
        self.holeNumberToBC = dict(zip(self.bestDS[:,0],
                                       self.bestDS[:,1:self.bestDS.shape[1]]))
             
        ## init the bcLabels according to the scoreMode. 
        if self.scoreMode == 'symmetric':
            bcLabels = [ makeBCLabel(a, b) for a,b in zip(self.barcodes, self.barcodes) ]
            bcLabels.append(NULL_BARCODE)
        elif self.scoreMode == 'asymmetric':
            bcLabels = [makeBCLabel(self.barcodes[i], self.barcodes[j]) 
                        for i in xrange(0, len(self.barcodes)) 
                        for j in xrange(i+1, len(self.barcodes))]
            bcLabels.append(NULL_BARCODE)
        elif self.scoreMode == 'paired':
            bcLabels = [ makeBCLabel(self.barcodes[i], self.barcodes[i+1]) 
                         for i in xrange(0, len(self.barcodes) - 1, 2) ]
            bcLabels.append(NULL_BARCODE)
        
        ## initialize the labels and the label map.
        self.bcLabels   = bcLabels
        self.labelToIdx = dict([ (v,k) for (k,v) in enumerate(self.bcLabels) ])

        ## it is a little confusing because, given a scoreMode then
        ## there is a new indexing system.
        self.scoreModeIdx = n.array([ self.getBarcodeTupleForZMW(hn)[0] 
                                      for hn in self.bestDS[:,0]])

        ## This will be used if we are using the extended "All" dataset.
        self.mapToAll = None

    def getBarcodeLabels(self):
        return self.bcLabels
    def getScoreMode(self):
        return self.scoreMode

    def nullBarcodeTuple(self):
        return (len(self.bcLabels) - 1, 0, 0)

    def getExtendedBarcodeInfoForZMW(self, holeNumber):
        # we'll do this on demand
        if not self.mapToAll:
            self.mapToAll = {}
            self.bcAllDS  = self.h5File[BC_DS_ALL_PATH][:]
            for r in xrange(0, self.bcAllDS.shape[0]):
                z = self.bcAllDS[r, 0]
                if not z in self.mapToAll:
                    self.mapToAll[z] = (r, r)
                else:
                    self.mapToAll[z] = (self.mapToAll[z][0], r)

        if holeNumber in self.mapToAll:
            extents = self.mapToAll[holeNumber]
            return self.bcAllDS[extents[0]:(extents[1]+1),:]
        else:
            return None


    def getBarcodeTupleForZMW(self, holeNumber):
        """Returns a tuple of (barcodeIdx, score, numberOfPasses)
        where barcodeIdx is an index into bcLabels"""
        try:
            ## XXX: something is not exactly right here w.r.t. the
            ## number of passes and the average scores.

            d = self.holeNumberToBC[holeNumber]
            if self.scoreMode == 'symmetric':
                return (d[1], d[2], d[0])

            elif self.scoreMode == 'asymmetric':
                (l, h) = (d[1], d[3]) if d[1] < d[3] else (d[3], d[1])
                idx = ((len(self.barcodes)-1)*l - l*(l-1)/2) + h-l-1
                return (idx, d[2] + d[4], d[0])

            elif self.scoreMode == 'paired':
                return (d[1]/2, d[2] + d[4], d[0])
        except:
            return self.nullBarcodeTuple()

    def getZMWsForBarcode(self, barcodeName):
        """Returns all the ZMWs that had barcodeName mapping to it,
        throws a BarcodeIdxException if there aren't any ZMWs for this
        barcode."""
        bcID = self.labelToIdx[barcodeName]
        msk = bcID == self.scoreModeIdx
       
        # if msk is all false then it won't return a 0-row data
        # structure. Furthermore, I can't return None because then
        # 'if' isn't happy because in some cases it is an array and in
        # others it is None. Have to raise an Exception.
        if n.any(msk):
            return self.bestDS[msk,:]
        else:
            raise BarcodeIdxException()
