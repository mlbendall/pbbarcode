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

#!/usr/bin/env python
import os
import sys
import argparse
import logging
import tempfile
import shutil
import pkg_resources

from pbcore.util.ToolRunner import PBMultiToolRunner

__version__ = ".01"

class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))

        subparsers = self.getSubParsers()

        desc = ['label-regions adds annotation to a region h5 file.']
        parser_m = subparsers.add_parser('label-regions', help = "Label regions from a region h5 file",
                                         description = "\n".join(desc), parents = [self.parser])
        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta',
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.h5',
                              help = 'Input region (rgn.h5), base (bas.h5), or fofn file containing region information')
        parser_m.add_argument('--prefix', metavar = 'file-name', dest = 'barcodeH5', default = None, 
                                 help = 'If prefix is specified then file-name.h5 is the resultant output')
        parser_m.add_argument('--overwrite', action='store_true', default = False,  
                              help = 'If overwrite is specified then the region information is updated in-place.')
        
        desc = ['label-reads adds information about barcode alignments to a cmp.h5 file',
                'from a previous call to "label-regions".']
        parser_s = subparsers.add_parser('label-reads', description = "\n".join(desc),
                                         parents = [self.parser], help = "Label reads from a barcode or region h5 file")
        parser_s.add_argument('inputFofn', metavar = 'barcode.fofn',
                              help = 'input barcode or region fofn file')
        parser_s.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                              help = 'cmp.h5 file to add barcode labels')
        
    def getVersion(self):
        return __version__
    
    def labelRegions(self):
        logging.debug("labeling regions.")
        if _isFOFN(self.args.inputFile):
            lines = open(self.args.inputFile, 'r').read().splitlines()
        else:
            lines = [ self.args.inputFile ]

        labels = [ BarcodeLabeler(self.args.barcodeFile, f).label() for line in lines ] 
        if (overwrite):
            # move barcode info into region file.
            pass
        else:
            pass
            # write barcode info to new barcode file.

    def labelReads(self):
        logging.debug("labeling reads.")
        

    def run(self):
        logging.debug("Arguments" + str(self.args))
        if self.args.subName == 'label-regions':
            self.labelRegions()
        elif self.args.subName == 'label-reads':
            self.labelReads()

if __name__ == '__main__':    
    sys.exit(Pbbarcode().start())
    


        

