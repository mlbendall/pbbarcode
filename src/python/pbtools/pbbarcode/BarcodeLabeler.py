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

from pbcore.io import FastaIO
from pbtools.pbbarcode import SWaligner

class RegionFile(object):
    HoleNumber = 0
    TypeIndex  = 1
    Start      = 2
    End        = 3
    Adapter    = 0
    Insert     = 1
    HQRegion   = 2

    def __init__(self, h5File):
        rgnTbl = h5File['PulseData/Regions']
        rgnDta = rgnTbl[:]
        rgnDta = rgnDta[np.lexsort((rgnDta[:,RegionFile.End], rgnDta[:,RegionFile.Start], 
                                    rgnDta[:,RegionFile.HoleNumber])),]
        zmwMap = {} 
        cZmw = rgnDta[0, RegionFile.HoleNumber]
        zmwMap[cZmw] = 0
        for i in xrange(1, rgnDta.shape[0]):
            if rgnDta[i, self.HoleNumber] != cZmw:
                newZmw = rgnDta[i, RegionFile.HoleNumber]
                oldZmw = rgnDta[i - 1, RegionFile.HoleNumber]
                zmwMap[newZmw] = i
                zmwMap[oldZmw] = (zmwMap[oldZmw], i)
                cZmw = newZmw
        zmwMap[cZmw] = (zmwMap[cZmw], rgnDta.shape[0])
        self.rgnDta = rgnDta
        self.zmwMap = zmwMap

    def zmws(self):
        for z in self.zmwMap.values():
            yield self.rgnDta[z[0]:z[1],]
           
class BasH5(object):
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

  
class BarcodeAligner(object):
    def __init__(self, barcodeFile, maxSequenceLength = 64):
        self.barcodes = dict([(x.getTag(), (x.sequence).upper()) for x in FastaIO.SimpleFastaReader(barcodeFile)])
        self.cMap     = dict(zip('ACGTacgt-N','TGCAtgca-N'))
        self.aligner  = SWaligner.SWaligner()
        
    def rc(self, s):
        return "".join([self.cMap[c] for c in s[::-1]])

    def align(self, seq1, seq2):
        return self.aligner.score(seq1,seq2)

    def getAlignmentScore(self, barcode, sequence):
        return max(self.align(barcode, sequence), self.align(self.rc(barcode), sequence))
    
    def scoreBarcodes(self, sequence):
        return [(name, self.getAlignmentScore(bcode, sequence)) for (name, bcode) in self.barcodes.items()]
        
class BarcodeLabeler(object):
    def __init__(self, barcodeFasta, basFile, rgnFile = None):
        bas = h5.File(basFile, 'r')
        if not rgnFile:
            rgn = bas
        else:
            rgn = h5.File(rgnFile, 'r')
        self.rgnFile = RegionFile(rgn)
        self.basFile = BasH5(bas)
        self.barcodes = BarcodeAligner(barcodeFasta)

    def processZMWs(self):
	pass        

    def processZMW(self, zmwBlk):
        w = ((zmwBlock[:,RegionFile.TypeIndex] == RegionFile.Adapter) | (zmwBlock[:,RegionFile.TypeIndex] == RegionFile.Insert))
        
        
        
if __name__ == "__main__":    
    bcl = BarcodeLabeler('/home/UNIXHOME/jbullard/projects/software/bioinformatics/tools/pbbarcode/etc/barcode.fasta',
                         '/mnt/data/vol18/2351565/0002/Results_SimpleP2B/m111121_220639_Spa_p1_b20.pls.h5')
    ziter = bcl.rgnFile.zmws()
    from IPython.Shell import IPShellEmbed; IPShellEmbed(argv=[])()

    
    
