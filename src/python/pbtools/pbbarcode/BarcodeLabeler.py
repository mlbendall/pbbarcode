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

from pbtools.pbbarcode.BarcodeH5Reader import LabeledZmw, \
    BARCODE_DELIMITER

__RC_MAP__ = dict(zip('ACGTacgt-N','TGCAtgca-N'))

class BarcodeScorer(object):
    def __init__(self, basH5, barcodeFasta,  
                 adapterSidePad = 0, insertSidePad = 4, 
                 scoreMode = 'symmetric', maxHits = 10, 
                 scoreFirst = False, startTimeCutoff = 1):
        """A BarcodeScorer object scores ZMWs and produces summaries
        of the scores. Various parameters control the behavior of the
        object, specifically the padding allows the user to add a
        little extra on each side of the adapter find for safety. The
        most relevant parameter is the scoreMode which dictates how
        the barcodes are scored, either paired or symmetric."""

        self.basH5 = basH5
        self.barcodeFasta = list(barcodeFasta)
        self.aligner = Aligner.SWaligner()
        self.barcodeLength = n.unique(map(lambda x : len(x.sequence), self.barcodeFasta))
        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

        self.barcodeSeqs = [(barcode.sequence.upper(), self._rc(barcode.sequence.upper())) 
                            for barcode in self.barcodeFasta]

        self.adapterSidePad = adapterSidePad
        self.insertSidePad = insertSidePad
        self.maxHits = maxHits

        if scoreMode not in ['symmetric', 'paired']:
            raise Exception("scoreMode must either be symmetric or paired")
        self._scoreMode = scoreMode

        self.scoreFirst = scoreFirst
        self.startTimeCutoff = startTimeCutoff

        self.forwardScorer = self.aligner.makeScorer([x[0] for x in self.barcodeSeqs])
        self.reverseScorer = self.aligner.makeScorer([x[1] for x in self.barcodeSeqs])
        
        logging.debug(("Constructed BarcodeScorer with scoreMode: %s, adapterSidePad: %d" + 
                       ", insertSidePad: %d, and scoreFirst: %r") % (scoreMode, adapterSidePad, 
                                                                     insertSidePad, scoreFirst))
        
    @property
    def movieName(self):
        return self.basH5.movieName
    
    def makeBCLabel(self, s1, s2):
        return BARCODE_DELIMITER.join((s1, s2))

    @property
    def barcodeLabels(self):
        """The barcode labels are function of the barcodeNames and the
        scoreMode, they represent the user-visible names."""
        if self.scoreMode == 'paired':
            return n.array([self.makeBCLabel(self.barcodeFasta[i].name,
                                             self.barcodeFasta[i+1].name) for i
                            in xrange(0, len(self.barcodeSeqs), 2)])
        else:
            return n.array([self.makeBCLabel(x.name, x.name) for x in self.barcodeFasta])

    @property
    def barcodeNames(self):
        """The barcode names are the FASTA names"""
        return n.array([x.name for x in self.barcodeFasta])

    @property
    def scoreMode(self):
        return self._scoreMode

    def _rc(self, s):
        return "".join([__RC_MAP__[c] for c in s[::-1]])
    
    def _flankingSeqs(self, zmw):
        def fromRange(rStart, rEnd):
            try:
                qSeqLeft = zmw.read(rStart - (self.barcodeLength + self.insertSidePad), 
                                    rStart + self.adapterSidePad).basecalls()
            except IndexError:
                qSeqLeft = None
            try:
                qSeqRight = zmw.read(rEnd - self.adapterSidePad, 
                                     rEnd + self.barcodeLength + 
                                     self.insertSidePad).basecalls()
            except IndexError:
                qSeqRight = None

            return (qSeqLeft, qSeqRight)

        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]
        
        seqs = [fromRange(start, end) for (start, end) in adapterRegions]

        # We only score the first barcode if we don't find any adapters
        # *and* the start time is less than the threshold. 
        scoredFirst = False
        if self.scoreFirst and not len(seqs):
            ## XXX: bug in pbcore which is throwing an exception here
            ## other than the one I'm catching.:
            try:
                s = zmw.zmwMetric('HQRegionStartTime')
                e = zmw.zmwMetric('HQRegionEndTime')
                # s<e => has HQ. 
                if s < e and s <= self.startTimeCutoff:
                    l = self.barcodeLength + self.insertSidePad
                    l = l if zmw.hqRegion[1] > l else zmw.hqRegion[1]
                    try:
                        bc = zmw.read(0, l).basecalls()
                        if len(bc) >= self.barcodeLength:
                            seqs.insert(0, (bc, None))
                            scoredFirst = True
                    except IndexError:
                        pass
            except:
                pass

        return (seqs, scoredFirst)

    def labelZmws(self, holeNumbers):
        """Return a list of LabeledZmws for input holeNumbers"""
        def scoreZmw(zmw):
            adapters, scoredFirst = self._flankingSeqs(zmw)
            adapterScores = [[]]*len(adapters)
            barcodeScores = n.zeros(len(self.barcodeSeqs))

            for i,adapter in enumerate(adapters):
                fscores  = self.forwardScorer(adapter[0])
                rscores  = self.reverseScorer(adapter[0])
                ffscores = self.forwardScorer(adapter[1])
                rrscores = self.reverseScorer(adapter[1])

                scored = 2 if adapter[0] and adapter[1] else \
                    1 if adapter[0] or adapter[1] else 0
                
                # An adapter score is the average barcode score for
                # each barcode -- that way, you can compare across
                # adapters even if the different adapters have
                # different numbers of flanking sequence. 
                if scored == 0:
                    adapterScores[i] = barcodeScores
                else:
                    adapterScores[i] = n.maximum((fscores + rrscores)/scored, 
                                                 (rscores + ffscores)/scored)

            barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
                else n.zeros(len(self.barcodeSeqs))

            return (zmw.holeNumber, len(adapters), barcodeScores, adapterScores,
                    scoredFirst)

        # o here is the record immediately above.
        def chooseSymmetric(o):
            p = n.argsort(-o[2])
            return LabeledZmw(o[0], o[1], p[0], o[2][p[0]], p[1], o[2][p[1]], o[3])
        def choosePaired(o):
            if o[1] == 1:
                s = n.array([max(o[2][i], o[2][i + 1]) for i in \
                                 xrange(0, len(self.barcodeSeqs), 2)])
                p = n.argsort(-s)
                s = s[p]
            else:
                # score the pairs by scoring the two alternate
                # ways they could have been put on the molecule. A
                # missed adapter will confuse this computation.
                scores  = o[3]
                results = n.zeros(len(self.barcodeSeqs)/2)
                for i in xrange(0, len(self.barcodeSeqs), 2):
                    pths = [0,0]
                    for j in xrange(0, len(scores)):
                        pths[j % 2] += scores[j][i]
                        pths[1 - j % 2] += scores[j][i + 1]
                        results[i/2] = max(pths)
                        
                p = n.argsort(-results)
                s = results[p]

            return LabeledZmw(o[0], o[1], p[0], s[0], p[1], s[1], o[3])
         
        if self.scoreMode == 'symmetric':
            choose = chooseSymmetric
        elif self.scoreMode == 'paired':
            choose = choosePaired
        else:
            raise Exception("Unsupported scoring mode in BarcodeLabeler.py")

        scored = [scoreZmw(self.basH5[zmw]) for zmw in holeNumbers] 
        return [choose(scoreTup) for scoreTup in scored if scoreTup[1]]


