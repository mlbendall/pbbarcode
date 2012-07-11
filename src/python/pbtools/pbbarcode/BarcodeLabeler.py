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
from pbcore.io.BasH5Reader import *
from pbcore.io.FastaIO import *
import pbtools.pbbarcode.SWaligner as Aligner
import numpy as n

_RC_MAP = dict(zip('ACGTacgt-N','TGCAtgca-N'))

class BarcodeScorer(object):
    def __init__(self, basH5, barcodeFasta, adapterSidePad = 0, insertSidePad = 4):
        self.basH5 = basH5
        self.barcodeFasta = list(barcodeFasta)
        self.aligner = Aligner.SWaligner()
        self.barcodeLength = n.unique(map(lambda x : len(x.sequence), self.barcodeFasta))
        self.barcodeSeqs = [(barcode.sequence.upper(), self.rc(barcode.sequence.upper())) 
                            for barcode in self.barcodeFasta]
        self.barcodeNames = n.array([barcode.name for barcode in self.barcodeFasta])
        self.adapterSidePad = adapterSidePad
        self.insertSidePad = insertSidePad
        
        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

    def score(self, s1, s2):
        if s1 and s2:
            return self.aligner.score(s1, s2)
        else:
            return 0
    
    def rc(self, s):
        return "".join([_RC_MAP[c] for c in s[::-1]])
    
    def flankingSeqs(self, zmw):
        def fromRange(rStart, rEnd):
            try:
                qSeqLeft = zmw.read(rStart - (self.barcodeLength + self.insertSidePad), 
                                    rStart + self.adapterSidePad).basecalls()
            except IndexError:
                qSeqLeft = None
            try:
                qSeqRight = zmw.read(rEnd - self.adapterSidePad, 
                                     rEnd + self.barcodeLength + self.insertSidePad).basecalls()
            except IndexError:
                qSeqRight = None
            return (qSeqLeft, qSeqRight)
        return [fromRange(start, end) for (start, end) in zmw.adapterRegions()]
        
    def scoreZMW(self, zmw):
        def scoreBarcode(barcode, adapterSeqs):
            return max(self.score(barcode[0], adapterSeqs[0]) + self.score(barcode[1], adapterSeqs[1]),
                       self.score(barcode[1], adapterSeqs[0]) + self.score(barcode[0], adapterSeqs[1]))

        barcodeScores = n.zeros(len(self.barcodeFasta))
        adapters = self.flankingSeqs(zmw)
        for adapter in adapters:
            barcodeScores += [scoreBarcode(barcode, adapter) for barcode in self.barcodeSeqs]
        
        return (zmw.holeNumber, len(adapters), barcodeScores)

    def scoreZMWs(self):
        return [self.scoreZMW(zmw) for zmw in self.basH5]

    def chooseBestBarcodes(self, s):
        def tabulate(o):
            p = n.argsort(-o[2])
            return (o[0], o[1], p[0], o[2][p[0]], p[1], o[2][p[1]])
        # if you have observed >= 1 adapter return something - this
        # still is lenient.
        return [tabulate(score) for score in s if score[1]] 
    

    
    
    
    
