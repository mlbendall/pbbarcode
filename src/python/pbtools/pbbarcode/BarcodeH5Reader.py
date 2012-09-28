import os
import logging

import h5py as h5
import numpy as n

from pbcore.io.BasH5Reader import *
from pbcore.io.CmpH5Reader import *

NULL_BARCODE   = "--"
BC_DS_PATH     = "BarcodeCalls/best"

def makeBCLabel(s1, s2):
    return NULL_BARCODE.join((s1, s2))

class BarcodeIdxException(Exception):
    pass

class BarcodeH5Reader(object):
    def __init__(self, fname):
        self.h5File    = h5.File(fname, 'r')
        self.bestDS    = self.h5File[BC_DS_PATH]
        self.scoreMode = self.bestDS.attrs['scoreMode']
        self.barcodes  = self.bestDS.attrs['barcodes']
        ## note: I hack of the ZMW.
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

    def getBarcodeTupleForZMW(self, holeNumber):
        """Returns a tuple of (barcodeIdx, score, numberOfPasses)
        where barcodeIdx is an index into bcLabels"""
        try:
            d = self.holeNumberToBC[holeNumber]
            if self.scoreMode == 'symmetric':
                return (d[1], d[2], d[0])
            elif self.scoreMode == 'asymmetric':
                (l, h) = (d[1], d[3]) if d[1] < d[3] else (d[3], d[1])
                idx = ((len(self.barcodes)-1)*l - l*(l-1)/2) + h-l-1
                return (idx, d[2]+d[4], d[0])
            elif self.scoreMode == 'paired':
                return (d[1]/2, d[2] + d[4], d[0])
        except:
            return (len(self.bcLabels) - 1, 0, 0)

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
