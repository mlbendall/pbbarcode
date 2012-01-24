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
from pbcore.io.cmph5 import factory
from pbtools.pbbarcode import SWaligner


__p4revision__ = "$Revision: #1 $"
__p4change__   = "$Change: 97803 $"
revNum         = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum      = int(__p4change__.strip("$").split(":")[-1])
__version__    = "%s-r%d-c%d" % ( pkg_resources.require("pbtools.barcode")[0].version, revNum, changeNum )
        
class BarcodeCallCopier(PBToolRunner):
    def __init__(self):
        desc = ['Command line tool for copying over barcode calls into cmp.h5 file.']
        super(BarcodeCallCopier, self).__init__('\n'.join(desc))
        self.parser.add_argument('barcodeFofn', metavar = 'barcode.fofn',
                                 help = 'input barcode fofn file')
        self.parser.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                                 help = 'cmp.h5 file to write to')

    def validateArgs(self):
        if not os.path.exists(self.args.barcodeFofn):
            self.parser.error('barcode.fofn file does file provided does not exist')
        if not os.path.exists(self.args.cmpH5):
            self.parser.error('aligned_reads.cmp.h5 file provided does not exist')
    
    def run(self):
        print "Input barcode.h5 file: %s" % self.args.barcodeFofn
        print "Input cmp.h5 file: %s" % self.args.cmpH5

        def makeKey(x):
            return os.path.basename(x).replace('.rgn.h5', '').replace('.bc.h5', '')
        def makeHash(x, name = 'Barcode'):
            h = h5.File(x)
            zmw = h['BarcodeData/Best/ZMW'][:]
            bcs = h['BarcodeData/Best/' + name][:]
            return dict(zip(zmw, bcs))
        
        barcodeFiles = open(self.args.barcodeFofn).read().splitlines()
        barcodeHash = dict([(makeKey(x), makeHash(x)) for x in barcodeFiles])
        scoreHash = dict([(makeKey(x), makeHash(x, "Score")) for x in barcodeFiles])

        cmph5 = factory.create(self.args.cmpH5, 'a')
        alnIndex = cmph5.asRecArray()
        movieDict = cmph5['/MovieInfo'].asDict('ID', 'Name')

        calls  = [""] * alnIndex.shape[0]
        scores = np.zeros(alnIndex.shape[0])

        for i in xrange(0, alnIndex.shape[0]):
            mid = alnIndex['MovieID'][i]
            zmw = alnIndex['HoleNumber'][i]
            try:
                calls[i]  = barcodeHash[movieDict[mid]][zmw]
                scores[i] = scoreHash[movieDict[mid]][zmw]
            except:
                calls[i] = ''

        def create(file, name, data, dtype):
            if name in file:
                del file[name]
            return file.create_dataset(name, data = data, dtype = dtype, maxshape = None)

        create(cmph5.h5File, "/AlnInfo/barcode", calls, dtype = h5.new_vlen(str))
        create(cmph5.h5File, "/AlnInfo/barcodeScore", scores, dtype = 'int32')
        
        # from IPython.Shell import IPShellEmbed; IPShellEmbed(argv=[])()

        
if __name__=="__main__":    
    sys.exit(BarcodeCallCopier().start())
    

      
        #     def makeKey(movie):
        #     return os.path.basename(y).replace('.pls.h5', '').replace('.bas.h5', '')
            
        # ## load the mappingCmpH5
        # cmph5 = factory.create(files.alignedReads.path)
        # alnIndex = cmph5.asRecArray()
        # movieDict = cmph5['/MovieInfo'].asDict('ID', 'Name')

        # ## bc dict
        # bcDict = {}
        # for movie in open(barcodeFofn).read().splitlines():
        #     h = H5.File(movie)
        #     zmw = h['BarcodeData/Best/ZMW'][:]
        #     bcs = h['BarcodeData/Best/Barcode'][:]
        #     bcDict[makeKey(movie)] = dict(zip(zmw, bcs))

        # print [ bcDict[movieDict[m]][h] for (m, h) in zip(alnIndex['MovieID'], alnIndex['HoleNumber']) ]
                
        # rgn = h5.File(self.args.rgnH5File, 'r')
        # bas = h5.File(self.args.basH5File, 'r')

        # self.rgnFile = RegionFile(rgn)
        # self.basFile = BasH5(bas)
        # self.barcodes = BarcodeAligner(self.args.barcodeFile)

        # res = [self.processZMW(x) for x in self.rgnFile.adapters()]
                
        # ## A semi-complete data structure. The maximum scoring barcode
        # ## for each adapter find.
        # scores = np.hstack([r[0] for r in res])
        # labels = np.hstack([r[1] for r in res])
        # zmws   = np.hstack([r[2] for r in res])

        # uZMWs     = np.sort(np.unique(zmws))
        # bestCalls = [""]*len(uZMWs)

        # def maxElt(v):
        #     u = np.unique(v)
        #     c = [0]*len(u)
        #     for i in xrange(0, len(u)):
        #         c[i] = sum(u[i] == v)
        #     return u[np.argmax(c)]

        # best = [(int(z), maxElt(labels[zmws == z])) for z in uZMWs] 
        # bestZmw = [b[0] for b in best]
        # bestBcode = [b[1] for b in best]
        
        # ## write the best data structure and the semi-complete and be done.
        # def create(file, name, data, dtype):
        #     if name in file:
        #         del file[name]
        #     return file.create_dataset(name, data = data, dtype = dtype, maxshape = None)

        # if self.args.barcodeH5:
        #     ofile = h5.File(self.args.barcodeH5)
        # else:
        #     ofile = rgn

        # create(ofile, "BarcodeData/Best/ZMW", data = bestZmw, dtype = "int32")
        # create(ofile, "BarcodeData/Best/Barcode", data = bestBcode, dtype = h5.new_vlen(str))
        # create(ofile, "BarcodeData/All/ZMW", data = zmws, dtype = "int32")
        # create(ofile, "BarcodeData/All/Score", data = scores, dtype = "int32")
        # create(ofile, "BarcodeData/All/Barcode", data = labels, dtype = h5.new_vlen(str))
        
        # try:
        #     ofile.close()
        #     rgn.close()
        #     bas.close()
        # except:
        #     pass
        

