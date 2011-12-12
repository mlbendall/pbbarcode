#!/usr/bin/env python
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
        self.expand = 20

        
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
        pass


    def run(self):
        print "input bas.h5 file: %s" % self.args.basH5File
        print "input rgn.h5 file: %s" % self.args.rgnH5File
        print "input barcode.fasta file: %s" % self.args.barcodeFile
        
        # 1.) make a copy of the region file for writing. 
        # 2.) iterate over ZMWs returning a new region chunk.
        # 3.) write the new file with appropriate metadata
        
        rgn = h5.File(self.args.rgnH5File)
        bas = h5.File(self.args.basH5File)
        rgnFile = RegionFile(rgn)
        basFile = BasH5(bas)

        from IPython.Shell import IPShellEmbed; IPShellEmbed(argv=[])()



if __name__=="__main__":    
    sys.exit(BarcodeLabeler().start())
