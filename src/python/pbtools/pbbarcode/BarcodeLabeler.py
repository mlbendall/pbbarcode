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
import logging

from pbcore.io.BasH5Reader import *
from pbcore.io.FastaIO import *
import pbtools.pbbarcode.SWaligner as Aligner
import numpy as n

_RC_MAP = dict(zip('ACGTacgt-N','TGCAtgca-N'))

class BarcodeScorer(object):
    def __init__(self, basH5, barcodeFasta, adapterSidePad = 0, insertSidePad = 4, 
                 scoreMode = 'symmetric', maxHits = 10, scoreFirst = False, 
                 startTimeCutoff = 1):
        self.basH5 = basH5
        self.barcodeFasta = list(barcodeFasta)
        self.aligner = Aligner.SWaligner()
        self.barcodeLength = n.unique(map(lambda x : len(x.sequence), self.barcodeFasta))
        self.barcodeSeqs = [(barcode.sequence.upper(), self.rc(barcode.sequence.upper())) 
                            for barcode in self.barcodeFasta]
        self.barcodeNames = n.array([barcode.name for barcode in self.barcodeFasta])
        self.adapterSidePad = adapterSidePad
        self.insertSidePad = insertSidePad
        self.scoreMode = scoreMode
        self.maxHits = maxHits
        self.scoreFirst = scoreFirst
        self.startTimeCutoff = startTimeCutoff

        self.forwardScorer = self.aligner.makeScorer([x[0] for x in self.barcodeSeqs])
        self.reverseScorer = self.aligner.makeScorer([x[1] for x in self.barcodeSeqs])

        logging.debug(("Constructed BarcodeScorer with scoreMode: %s, adapterSidePad: %d" + 
                       ", and insertSidePad: %d") % (scoreMode, adapterSidePad, insertSidePad))
        
        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

    # different path for scoring - slower.
    def score(self, s1, s2):
        if s1 and s2:
            return self.aligner.score(s1, s2)
        else:
            return 0
    def _scoreZMW_v1(self, zmw):
        def scoreBarcode(barcode, adapterSeqs):
            return max(self.score(barcode[0], adapterSeqs[0]) + self.score(barcode[1], adapterSeqs[1]),
                       self.score(barcode[1], adapterSeqs[0]) + self.score(barcode[0], adapterSeqs[1]))
    
        barcodeScores = n.zeros(len(self.barcodeFasta))
        adapters = self.flankingSeqs(zmw)
        for adapter in adapters:
            barcodeScores += [scoreBarcode(barcode, adapter) for barcode in self.barcodeSeqs]
        
        return (zmw.holeNumber, len(adapters), barcodeScores)
    
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

        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]
        
        seqs = [fromRange(start, end) for (start, end) in adapterRegions]

        ## try to score the first barcode.
        if self.scoreFirst:
            s = zmw.zmwMetric('HQRegionStartTime')
            e = zmw.zmwMetric('HQRegionEndTime')
            # s<e => has HQ. 
            if s < e and s <= self.startTimeCutoff:
                l = self.barcodeLength + self.insertSidePad
                l = l if zmw.hqRegion[1] > l else zmw.hqRegion[1]
                seqs.insert(0, (zmw.read(0, l).basecalls(), None))
        return seqs

    def scoreZMW(self, zmw):
        adapters = self.flankingSeqs(zmw)
        adapterScores = [[]]*len(adapters)
        barcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            fscores = self.forwardScorer(adapter[0])
            rscores = self.reverseScorer(adapter[0])
            ffscores = self.forwardScorer(adapter[1])
            rrscores = self.reverseScorer(adapter[1])
            adapterScores[i] = n.maximum(fscores + rrscores, 
                                         rscores + ffscores)
            
        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else n.zeros(len(self.barcodeSeqs))
        
        return (zmw.holeNumber, len(adapters), barcodeScores, adapterScores)

    def scoreZMWs(self, zmws = None):
        if zmws is None:
            zmws = self.basH5.sequencingZmws
        
        return [self.scoreZMW(self.basH5[zmw]) for zmw in zmws]

    # The expected record that is returned is: 
    # 
    # (holeNumber, nAdapters, barcodeIdx1, barcodeScore1, barcodeIdx2, barcodeScore2)
    # 
    # XXX : I'm not using named tuples because of performance, but readability truly suffers. 
    def chooseBestBarcodes(self, allScores):
        # remove cases where you obvserve 0 adapters.
        fScores = filter(lambda x : x[1], allScores)

        if self.scoreMode in ['symmetric', 'asymmetric']:
            def tabulate(o):
                p = n.argsort(-o[2])
                return (o[0], o[1], p[0], o[2][p[0]], p[1], o[2][p[1]])
        elif self.scoreMode in ['paired']:
            def tabulate(o):
                # bug 22035 - the problem is that I score both F1,R1
                # on both sides of the molecule. In the case when you
                # only see one adapter, then you are going to sum both
                # R1+F1 on the same putative barcode - this might be a
                # lower score than some alternate F or
                # R.
                if o[1] == 1:
                    p = n.argsort(-o[2])[0]
                    # now p can be either F1 or R1, but below, I need
                    # it to be F1 so if it is odd, subtract 1. And we
                    # need to mod it by 2 because our selection below
                    # will multiply it up.
                    p = (p if p % 2 == 0 else p - 1)/2
                else:
                    # score the pairs by scoring the two alternate
                    # ways they could have been put on the molecule. A
                    # missed adapter is an issue. 
                    scores  = o[3]
                    results = n.zeros(len(self.barcodeSeqs))

                    for i in xrange(0, len(self.barcodeSeqs), 2):
                        pths = [0,0]
                        for j in xrange(0, len(scores)):
                            pths[j % 2] += scores[j][i]
                            pths[1 - j % 2] += scores[j][i + 1]

                        results[i] = max(pths)
                    p = n.argmax(results)/2

                    # This block was how I used to score the pairs.
                    # p = n.argsort([-(o[2][i] + o[2][i+1]) for i in 
                    #                 xrange(0, len(o[2]) - 1, 2)])[0]

                return (o[0], o[1], 2*p, o[2][2*p], 2*p + 1, o[2][2*p + 1])

        else:
            raise Exception("Unsupported scoring mode in BarcodeLabeler.py")

        return [tabulate(score) for score in fScores]
