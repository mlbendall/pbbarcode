import sys
import os
import shutil
import glob
import logging
import tempfile
import subprocess
import pkg_resources
import h5py as h5
import numpy as np

from pbcore.util.ToolRunner import PBToolRunner
from pbcore.io import FastaIO

__p4revision__ = "$Revision: #1 $"
__p4change__   = "$Change: 97803 $"
revNum         = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum      = int(__p4change__.strip("$").split(":")[-1])
__version__    = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.barcode")[0].version, revNum, changeNum )

class RegionFile(object):
    def __init__(self, h5File):
        ## XXX: don't need to sort, you can easily make lookup constant.
        ## XXX: need to read the following indices from the attributes.
        self.TypeIndex  = 1
        self.HoleNumber = 0
        self.Start      = 2
        self.End        = 3
        self.Adapter    = 0
        rgnTbl = h5File['PulseData/Regions']
        rgnDta = rgnTbl[:]
        rgnDta = rgnDta[rgnDta[:,self.TypeIndex] == self.Adapter,]
        self.rgnDta = rgnDta[np.lexsort((rgnDta[:,self.End], rgnDta[:,self.Start], rgnDta[:,self.HoleNumber])),]

    def getAdaptersForHole(self, holeNumber):
        return self.rgnDta[self.rgnDta[:,self.HoleNumber] == holeNumber,:]
    
    def adapters(self):
        for z in np.unique(self.rgnDta[:,self.HoleNumber]):
            yield self.getAdaptersForHole(z)
            
class BasH5(object):
    ## XXX: this has to be built on basH5 class and work with both CCS and raw reads.
    def __init__(self, h5File, baseCallRoot = '/PulseData/BaseCalls'):
        self.baseCallsDG = h5File[baseCallRoot]
        self.baseCalls = self.baseCallsDG['Basecall']
        self.ZMWDG = self.baseCallsDG['ZMW']
        self.numEvent = self.ZMWDG['NumEvent'][:]
        self.holeNumber = self.ZMWDG['HoleNumber'][:]
        self.holeStatus = self.ZMWDG['HoleStatus'][:]
        self.hn2status = dict(zip(self.holeNumber, self.holeStatus))
        holeType = self.ZMWDG['HoleStatus'].attrs['LookupTable']
        self._holeTypeDict = dict(zip(range(len(holeType)), list(holeType)))

        endOffset = np.cumsum(self.numEvent)
        beginOffset = np.hstack( (np.array([0]), endOffset[0:-1] ) )
        offsets = np.vstack ( (beginOffset, endOffset) )

        self.hn2offsets = dict(zip(self.holeNumber, zip(*offsets)))
        self.baseCallsDG['Basecall']

    def getSequence(self, hole, start, end):
        zstart, zend = self.hn2offsets[hole]
        return self.baseCalls[(zstart + start):(zstart + end)]

    def getSequenceAsString(self, hole, start, end):
        return ''.join([chr(c) for c in self.getSequence(hole, start, end)])


def align(seq1, seq2, penalty = -1, match = 2):
    smat = np.zeros((len(seq1)+1, len(seq2)+1))
    bestScore = 0
    for i in xrange(1, smat.shape[0]):
        for j in xrange(1, smat.shape[1]):
            iscore = smat[i,j-1] + penalty
            dscore = smat[i-1,j] + penalty
            mscore = smat[i-1,j-1] + (match if seq1[i-1] == seq2[j-1] else penalty)
            smat[i,j] = np.max((0, iscore, dscore, mscore))
            if smat[i,j] >= bestScore:
                bestScore = smat[i,j]
    return bestScore

cMap = dict(zip('ACGTacgt-N','TGCAtgca-N'))
def rc(s):
    return "".join([cMap[c] for c in s[::-1]])
    
class BarcodeAligner(object):
    def __init__(self, barcodeFile, maxSequenceLength = 64):
        self.barcodes = dict([(x.getTag(), (x.sequence).upper()) for x in FastaIO.SimpleFastaReader(barcodeFile)])

    def getAlignmentScore(self, barcode, sequence):
        ## XXX: this bit will have to be fast. Options include having pre-allocated matrices for each barcode.
        return max(align(barcode, sequence), align(rc(barcode), sequence))

        
    def scoreBarcodes(self, sequence):
        return [(name,self.getAlignmentScore(bcode, sequence)) for (name, bcode) in self.barcodes.items()]
        

class BarcodeLabeler(PBToolRunner):
    def __init__(self):
        desc = ['Tool for labeling ZMWs with barcodes and creating a new regions file excluding the regions for future analysis']
        super(BarcodeLabeler, self).__init__('\n'.join(desc))
        self.parser.add_argument('basH5File', metavar='input.bas.h5',
                                 help='input .bas.h5 filename')
        self.parser.add_argument('rgnH5File', metavar='input.rgn.h5',
                                 help='input .rgn.h5 filename')
        self.parser.add_argument('barcodeFile', metavar='barcode.fasta',
                                 help='input barcode fasta file')
        self.expand = 30

        
    ## XXX: I don't think I should have to override this.
    def getVersion(self):
        return __version__

    def validateArgs(self):
        if not os.path.exists(self.args.basH5File):
            self.parser.error('input.bas.h5 file provided does not exist')
        if not os.path.exists(self.args.rgnH5File):
            self.parser.error('input.rgn.h5 file provided does not exist')
        if not os.path.exists(self.args.barcodeFile):
            self.parser.error('barcode.fasta file provided does not exist')
    
    def processZMW(self, npBlock):
        ## this method returns a new ZMW block including adapter hits.
        scores = np.zeros(shape = npBlock.shape[0])
        labels = [""]*npBlock.shape[0]
        ZMWs = np.zeros(shape = npBlock.shape[0])
        
        for j in xrange(0, npBlock.shape[0]):
            zmw,start,end = npBlock[j, [0,2,3]]
            
            ## score the lh bc.
            rstart = start - self.expand if start > self.expand else 0
            bc1 = self.barcodes.scoreBarcodes(self.basFile.getSequenceAsString(zmw,rstart,start))

            ## score the rh bc.
            rend = end + self.expand ## XXX: the bound on this is not obvious.
            bc2 = self.barcodes.scoreBarcodes(self.basFile.getSequenceAsString(zmw,end,rend))
        
            ## sum the two scores.
            bcscores = [(y[0][0], y[0][1] + y[1][1]) for y in zip(bc1,bc2)]
            labels[j],scores[j] = bcscores[np.argmax([x[1] for x in bcscores])]
            ZMWs[j] = zmw
        return (scores, labels, ZMWs)

    ##
    ## Things to think about:
    ## 1.) Should we "shrink" the insert region corresponding to the
    ##     barcode size - if we do this, then asymetric designs, i.e.,
    ##     adapter-bc-ins-adpater will have problems.
    ##
    ## 2.) Should we return some "different" datastructure?
    ##
    def run(self):
        print "Input bas.h5 file: %s" % self.args.basH5File
        print "Input rgn.h5 file: %s" % self.args.rgnH5File
        print "Input barcode.fasta file: %s" % self.args.barcodeFile
 
        rgn = h5.File(self.args.rgnH5File)
        bas = h5.File(self.args.basH5File)

        self.rgnFile = RegionFile(rgn)
        self.basFile = BasH5(bas)
        self.barcodes = BarcodeAligner(self.args.barcodeFile)

        from IPython.Shell import IPShellEmbed; IPShellEmbed(argv=[])()
        print "Done!"

        for x in self.rgnFile.adapters():
            print self.processZMW(x)

            # res = [self.processZMW(x) for x in self.rgnFile.adapters()]
        
        ## A semi-complete data structure. The maximum scoring barcode
        ## for each adapter find.
        scores = np.hstack([r[0] for r in res])
        labels = np.hstack([r[1] for r in res])
        zmws   = np.hstack([r[2] for r in res])
                

if __name__=="__main__":    
    sys.exit(BarcodeLabeler().start())
