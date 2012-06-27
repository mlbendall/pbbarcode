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
from pbtools.pbbarcode import SWaligner

class BarcodeFofnWrapper(PBToolRunner):
    def __init__(self):
        desc = ['Wrapper to work with fofns']
        super(BarcodeFofnWrapper, self).__init__('\n'.join(desc))
        self.parser.add_argument('--out-fofn', metavar = 'barcode.fofn', dest = 'outFofn',
                                 default = '', help = 'output barcode fofn file')
        self.parser.add_argument('--in-fofn', metavar = 'input.fofn', dest = 'inFofn',
                                 default = '', help = 'input pls fofn file')
        self.parser.add_argument('--bc-fasta', metavar = 'bc.fasta', dest = 'bcFasta',
                                 default = '', help = 'barcode fasta file')
 
    def run(self):
        outFofn = open(self.args.outFofn, 'w')
        inFofn  = open(self.args.inFofn, 'r').read().splitlines()
        for pls in inFofn:
            bcName = os.path.basename(pls).replace('.pls.h5', '.bc.h5').replace('.bas.h5', '.bc.h5')
            bcFile = os.path.dirname(self.args.outFofn) + '/' + bcName
            os.system("barcodeLabeler.py --out-file %s %s %s %s" % (bcFile, pls, pls, self.args.bcFasta))
            outFofn.write(bcFile + "\n")

        outFofn.close()

    
    def validateArgs(self):
        if not os.path.exists(self.args.inFofn):
            self.parser.error('input.fofn does not exist')
        if not os.path.exists(self.args.bcFasta):
            self.parser.error('barcode.fasta file provided does not exist')
            
if __name__=="__main__":    
    sys.exit(BarcodeFofnWrapper().start())


