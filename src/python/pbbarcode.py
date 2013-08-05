#!/usr/bin/env python
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
import os
import sys
import argparse
import logging
import tempfile
import shutil
import pkg_resources
import re
import subprocess
import random
import shutil

from multiprocessing import Pool

import h5py as h5
import numpy as n

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io.BasH5Reader import *
from pbcore.io.CmpH5Reader import *
from pbcore.io.BarcodeH5Reader import *
from pbcore.io import FastaReader, FastqWriter, FastqRecord, \
    FastaWriter, FastaRecord

from pbtools.pbbarcode.BarcodeLabeler import *
from pbtools.pbh5tools.CmpH5Utils import copyAttributes

__version__ = ".4"

# Paths to the Barcode Datasets in the cmp.h5 file.
BC_ALN_INFO_DS = "AlnInfo/Barcode"
BC_INFO_NAME   = "BarcodeInfo/Name"
BC_INFO_ID     = "BarcodeInfo/ID"

SCORE_MODES    = ['symmetric', 'paired']

BAS_PLS_REGEX = r'\.ba[x|s]\.h5$|\.pl[x|s]\.h5$|\.cc[x|s]\.h5$'
BARCODE_EXT   = '.bc.h5'
BC_REGEX      = r'\.bc\.h5'

def movieNameFromFile(fn):
    return re.sub('|'.join((BC_REGEX, BAS_PLS_REGEX)) , '', 
                  os.path.basename(fn))

def makeBarcodeH5FromBasH5(basH5):
    """The workhorse function for creating a barcode H5 file from a
    base H5 file."""
    labeler = BarcodeScorer(basH5, FastaReader(runner.args.barcodeFile),
                            runner.args.adapterSidePad, runner.args.insertSidePad,
                            scoreMode = runner.args.scoreMode, 
                            maxHits = runner.args.maxAdapters,
                            scoreFirst = runner.args.scoreFirst, 
                            startTimeCutoff = runner.args.startTimeCutoff)
    if runner.args.nZmws < 0:
        zmws = basH5.sequencingZmws
    else:
        zmws = basH5.sequencingZmws[0:runner.args.nZmws]

    logging.debug("Labeling %d ZMWs from: %s" % (len(zmws), basH5.filename))
    labeledZmws = labeler.labelZmws(zmws)
    logging.debug("Labeled %d ZMWs" % len(labeledZmws))
    
    outBase = re.sub(BAS_PLS_REGEX, BARCODE_EXT, 
                     os.path.basename(basH5.filename))
    outFile = '/'.join((runner.args.outDir, outBase))
    logging.debug("Writing to: %s" % outFile)

    writeBarcodeH5(labeledZmws, labeler, outFile, 
                   runner.args.saveExtendedInfo)
    return outFile

def mpWrapper(f):
    return makeBarcodeH5FromBasH5(BasH5Reader(f))

def makeBarcodeFofnFromBasFofn():
    inputFofn = runner.args.inputFile
    inFiles = open(inputFofn).read().splitlines()
    
    if not all(map(os.path.exists, inFiles)):
        raise IOError("All files in input.fofn must exist.")
    
    logging.debug("Using %d processes." % runner.args.nProcs)
    if runner.args.nProcs <= 1:
        newFiles = map(mpWrapper, inFiles)
    else:
        pool = Pool(runner.args.nProcs)
        newFiles = pool.map(mpWrapper, inFiles)

    oFile = open(runner.args.outFofn, 'w')
    for nF in newFiles:
        oFile.write(nF + "\n")
    oFile.close()

def labelAlignments():
    logging.info("Labeling alignments using: %s" % runner.args.inputFofn)
    bcFofn = BarcodeH5Fofn(runner.args.inputFofn)
    
    with CmpH5Reader(runner.args.cmpH5) as cmpH5:
        bcDS = n.zeros((len(cmpH5), 5), dtype = "int32")

        for (i, aln) in enumerate(cmpH5):
            bcReader = bcFofn.readerForMovie(aln.movieInfo.Name)
            try:
                lZmw = bcReader.labeledZmwFromHoleNumber(aln.HoleNumber)
                if lZmw.nScored < runner.args.minNumBarcodes or \
                        lZmw.averageScore < runner.args.minAvgBarcodeScore or \
                        lZmw.scoreRatio < runner.args.minScoreRatio:
                    lZmw = None
            except KeyError:
                lZmw = None
            
            if lZmw:
                bcDS[i,:] = n.array([lZmw.nScored, lZmw.bestIdx, lZmw.bestScore,
                                     lZmw.secondBestIdx, lZmw.secondBestScore])
            else:
                # either no barcode was found for this guy or they got
                # filtered, hence the NULL_BARCODE
                bcDS[i,:] = n.array([0, 
                                     len(bcReader.barcodeLabels), 0, 
                                     len(bcReader.barcodeLabels), 0])

    # write to the cmp.h5 file.
    H5 = h5.File(runner.args.cmpH5, 'r+')
    if BC_INFO_ID in H5:
        del H5[BC_INFO_ID]
    if BC_INFO_NAME in H5:
        del H5[BC_INFO_NAME]

    # we use the first one to get the labels, if somehow they
    # don't have all of the same stuff that will be an issue.
    bcLabels = n.concatenate((bcFofn.barcodeLabels, n.array([BARCODE_DELIMITER]))) 
    H5.create_dataset(BC_INFO_ID, data = n.array(range(0, len(bcLabels))), 
                      dtype = 'int32')
    H5.create_dataset(BC_INFO_NAME, data = bcLabels, dtype = h5.new_vlen(str))
    if BC_ALN_INFO_DS in H5:
        del H5[BC_ALN_INFO_DS]
    bcDS = H5.create_dataset(BC_ALN_INFO_DS, data = bcDS, dtype = 'int32')
    bcDS.attrs['ColumnNames'] = n.array(['count', 'index1', 'score1', 'index2', 
                                         'score2'])
    bcDS.attrs['BarcodeMode'] = bcFofn.scoreMode
    H5.close()

def zipFofns(*inFofns):
    """Take inputFofns and return n tuples of length len(inFofns)
    where n is the number of entries in each FOFN."""
    def readAndSort(inFile):
        lines = n.array(open(inFile).read().splitlines())
        lines = lines[n.array(n.argsort([movieNameFromFile(fofnLine) for 
                                         fofnLine in lines]))]
        return lines

    sortedFofns = [readAndSort(inFofn) for inFofn in inFofns]
    l = map(len, sortedFofns)
    if len(n.unique(l)) != 1:
        raise Exception("Fofns don't match, unequal number of inputs.")
    else:
        for i in xrange(0, n.unique(l)):
            if len(n.unique([movieNameFromFile(sortedFofn[i]) for 
                             sortedFofn in sortedFofns])) != 1:
                raise Exception("Fofn elements don't match, movies differ.")
    
    # need to un-arrayify these guys
    return zip(*map(list, sortedFofns))

def filterZmws(zmwsForBCs):
    """Apply various filterings passed by the user. There are somewhat
    different semantics for CCS filtering and subread filtering in
    terms of the raw primary metrics available, e.g.,
    HQRegionStartTime is unavailable for the CCS data and somewhat
    irrelevant."""
    def getHQStart(zmw):
        try:
            return zmw.zmwMetric('HQRegionStartTime')
        except:
            return 0
        
    def getReadScore(zmw):
        return zmw.zmwMetric("ReadScore")

    def molLenGuess(zmw):
        if zmw.baxH5.hasRawBasecalls:
            return max(map(len, zmw.subreads)) if zmw.subreads else 0
        else:
            return len(zmw.ccsRead) if zmw.ccsRead else 0

    def zmwFilterFx(tup):
        zmw, lZmw = tup

        mlGuess = molLenGuess(zmw)
        if not mlGuess:
            return False

        avgScore    = lZmw.averageScore
        numScored   = lZmw.nScored
        scoreRatio  = lZmw.scoreRatio
        hqStart     = getHQStart(zmw)
        readScore   = getReadScore(zmw)

        ## XXX : still need to detect the chimeras
        if mlGuess < runner.args.minMaxInsertLength or \
                hqStart > runner.args.hqStartTime or \
                readScore < runner.args.minReadScore or \
                avgScore < runner.args.minAvgBarcodeScore or \
                numScored < runner.args.minNumBarcodes or \
                scoreRatio < runner.args.minScoreRatio:
            return False
        else:
            return True

    return { k:filter(zmwFilterFx, v) for k,v in zmwsForBCs.items() }

def _warnOnce():
    var = []
    def warnOnce(msg):
        if not var:
            logging.warn(msg)
        var.append(1)
    return warnOnce
warnOnce = _warnOnce()

def getFastqRecords(zmw, lZmw = None):
    if zmw.baxH5.hasRawBasecalls and zmw.baxH5.hasConsensusBasecalls:
        # Only examine this parameter when passed both.
        if runner.args.subreads:
            reads = zmw.subreads
        else:
            reads = [zmw.ccsRead]
    elif zmw.baxH5.hasRawBasecalls:
        if runner.args.subreads:
            warnOnce("`subreads` argument is ignored when using >= 2.1" + 
                     "bas.h5 data as input.")
        reads = zmw.subreads
    else:
        if runner.args.subreads:
            warnOnce("`subreads` argument is ignored when using >= 2.1" + 
                     "ccs.h5 data as input.")
        reads = [zmw.ccsRead]

    extra = (" %g %g" % (round(zmw.zmwMetric("ReadScore"), 2), 
                         round(lZmw.averageScore, 2))) if lZmw else ""
    
    return [FastqRecord(read.readName + extra,
                        read.basecalls(),
                        read.QualityValue()) for read in reads if read]

def getFastqs():
    zmwsByBarcode = getZmwsForBarcodes()
    logging.debug("Pre-filter: Average number of ZMWs per barcode: %d" % 
                  n.mean([len(zmwsByBarcode[k]) for k in zmwsByBarcode.keys()]))

    zmwsByBarcode = filterZmws(zmwsByBarcode) 
    logging.debug("Post-filter: Average number of ZMWs per barcode: %d" % 
                  n.mean([len(zmwsByBarcode[k]) for k in zmwsByBarcode.keys()]))
    
    def getReadData(zmws):
        recs = [getFastqRecords(zmw,lZmw) for zmw,lZmw in zmws]
        recs = filter(lambda x : x, recs)
        return [elt for sublst in recs for elt in sublst]

    return {k:getReadData(zmws) for k, zmws in zmwsByBarcode.iteritems()}

def emitFastqs():
    outFiles = getFastqs()
    outDir   = runner.args.outDir
    fasta    = runner.args.fasta

    if runner.args.unlabeledZmws:
        outFiles['UNLABLED'] = getUnlabeledZmws()

    if not os.path.exists(runner.args.outDir):
        os.makedirs(runner.args.outDir)

    if fasta:
        writer = FastaWriter
        def record(n, s, qv):
            return FastaRecord(n, s)
    else:
        writer = FastqWriter
        record = FastqRecord
    
    l = 'a' if runner.args.fasta else 'q'
    for k in outFiles.keys():
        if outFiles[k]:
            with writer("%s/%s.fast%s" % (runner.args.outDir, k, l)) as w:
                for e in outFiles[k]:
                    tlen = len(e.sequence)-runner.args.trim
                    r = record(e.name, e.sequence[runner.args.trim:tlen],
                               e.quality[runner.args.trim:tlen])
                    if r:
                        w.writeRecord(r)

def getUnlabeledZmws():
    """Return FASTQ records for ZMWs which do not have a barcode label"""
    unlabeledZmws = []

    for basFile, barcodeFile in zipFofns(runner.args.inputFofn, 
                                         runner.args.barcodeFofn):
        basH5 = BasH5Reader(basFile)
        bcH5  = BarcodeH5Reader(barcodeFile)
        sdiff = basH5.sequencingZmws[~n.in1d(basH5.sequencingZmws,  
                                             bcH5.labeledZmws.keys())]
        for hn in sdiff:
            unlabeledZmws.append(basH5[hn])

    return reduce(lambda x,y : x+y, [getFastqRecords(unlabeledZmw) for 
                                     unlabeledZmw in unlabeledZmws])

def getZmwsForBarcodes():
    """dictionary of pbcore.io.Zmw and LabeledZmw indexed by barcode
    label"""
    zmwsForBCs = {} 

    for basFile, barcodeFile in zipFofns(runner.args.inputFofn, 
                                         runner.args.barcodeFofn):
        basH5 = BasH5Reader(basFile)
        bcH5  = BarcodeH5Reader(barcodeFile)
        for label in bcH5.barcodeLabels:
            lZmws = bcH5.labeledZmwsFromBarcodeLabel(label)
            for lZmw in lZmws:
                zmw = basH5[lZmw.holeNumber]
                if not label in zmwsForBCs.keys():
                    zmwsForBCs[label] = []
                zmwsForBCs[label].append((zmw, lZmw))

    return zmwsForBCs

def gconFunc(tp):
    # called bcause multiprocess
    rootDir, barcode = tp
    bcdir = "/".join((rootDir, barcode))

    ## call gcon
    logging.info("In gconFunc for: %s" % barcode)

    cmd = "gcon.py r --min_cov 3 %s/subreads.fasta %s/seed_read.fasta -d %s" % \
        (bcdir, bcdir, bcdir)
    subprocess.call(cmd, shell = True)

    ## check to see if the file is empty
    r = FastaReader("%s/g_consensus.fa" % bcdir)
    
    if not list(r)[0].sequence:
        return None
    
    ## setup the blasr / sam / quiver stuff.
    logging.info("Setup regions file, now running blasr through quiver.")

    cmd = ('blasr %s %s/g_consensus.fa -nproc 1 -sam -regionTable %s/region.fofn -out ' + \
               '%s/aligned_reads.sam') % (runner.args.inputFofn, bcdir, bcdir, bcdir)
    subprocess.call(cmd, shell = True)
    
    logging.info("Done with blasr in: %s" % barcode)
    
    cmd = 'samtoh5 %s/aligned_reads.sam %s/g_consensus.fa %s/aligned_reads.cmp.h5' % \
        (bcdir, bcdir, bcdir)
    subprocess.call(cmd, shell = True)
         
    logging.info("Done with samtoh5 in: %s" % barcode)
    
    cmd = ('loadPulses %s %s/aligned_reads.cmp.h5 -byread -metrics ' + \
               'QualityValue,InsertionQV,MergeQV,DeletionQV,DeletionTag,SubstitutionTag,' + \
               'SubstitutionQV') % (runner.args.inputFofn, bcdir)
    subprocess.call(cmd, shell = True)
    logging.info("Done with loadPulses in: %s" % barcode)
    
    cmd = 'cmph5tools.py sort --inPlace %s/aligned_reads.cmp.h5' % bcdir
    subprocess.call(cmd, shell = True)
    logging.info("Done with cmph5tools.py: %s" % barcode)

    cmd = ('quiver %s/aligned_reads.cmp.h5 --outputFilename %s/q_consensus.fasta ' + \
               '--referenceFilename %s/g_consensus.fa') % (bcdir, bcdir, bcdir)
    subprocess.call(cmd, shell = True)
    logging.info("Done with quiver: %s" % barcode)
	
        
    ## append results to output file.
    bcCons = "%s/%s/q_consensus.fasta" % (rootDir, barcode)
    if os.path.exists(bcCons):
        return FastaRecord(barcode, list(FastaReader(bcCons))[0].sequence)
    else:
        return None

def subsampleReads(e):
    logging.info("starting with %d zmws" % len(e))
    if runner.args.nZmws > 0:
        k = runner.args.nZmws if runner.args.nZmws < len(e) else len(e)    
    elif runner.args.subsample < 1:
        k = int(len(e)*runner.args.subsample)
    else:
        k = len(e)
    i = n.array(random.sample(range(0, len(e)), k), dtype = int)
    logging.info("subsampled down to: %d" % len(i))
    return [e[j] for j in i]

def callConsensus():
    def makeReadAndReads(zmwsForBC):
        ccsData = filter(lambda x:x, [zmw.ccsRead for _,_,zmw in zmwsForBC if zmw])
        srData  = reduce(lambda x,y : x+y, [zmw.subreads for zmw,_,_ in
                                            zmwsForBC if zmw], [])
        if not srData and not ccsData:
            return (None,None)

        def getSeedRead(reads, lq = 80, uq = 90):
            lens = map(len, reads)
            candidateRange = (n.percentile(lens, 80), 
                              n.percentile(lens, 90))
            pfReads = [read for read,l in zip(reads, lens) if 
                       l >= candidateRange[0] and l <= candidateRange[1]]
            pfReads.sort(key = lambda x : -x.zmw.readScore)
            return pfReads[0] if len(pfReads) else None

        if ccsData:
            # it is much safer to pick basically the longest CCS read
            # than the subread.
            seedRead = getSeedRead(ccsData, lq = 975, uq = 1)
            if not seedRead:
                seedRead = getSeedRead(srData)
                logging.info("Unable to use a CCS read for seed read.")
            else:
                logging.info("Using CCS read for seed read.")

            ccsZmws  = [x.holeNumber for x in ccsData]
            allReads = ccsData
            for sr in srData:
                if sr.holeNumber not in ccsZmws:
                    allReads.append(sr)
        else:
            logging.info("Using raw read for seed read")
            seedRead = getSeedRead(srData)
            allReads = srData
        
        return (seedRead, allReads)
    

    # check to make sure that you have the necessary dependencies,
    # i.e., hgap script, blasr, etc.
    try:
        import pbtools.pbdagcon
    except ImportError:
        raise ImportError("Unable to find dependency `pbdagcon` - please install.")

    # retrieve ZMWs by barcode
    zmwsForBCs = getZmwsForBarcodes()
    
    # subsample
    zmwsForBCs = {k:subsampleReads(v) for k,v in zmwsForBCs.items()}

    logging.info("unfiltered average zmws per barcode: %g" % 
                 n.round(n.mean(map(len, zmwsForBCs.values()))))

    # filter ZMWs
    zmwsForBCs = filterZmws(zmwsForBCs)
    
    logging.info("filtered average zmws per barcode: %g" % 
                 n.round(n.mean(map(len, zmwsForBCs.values()))))

    # now choose the best subread to seed the assembly
    if runner.args.ccsFofn:
        # XXX: This part depends on the filenames of the ccs and input
        # fofns, this is essentially a workaround to the fact the the
        # part isn't part of the API
        ccsReaders = {movieNameFromFile(l):BasH5Reader(l) for l in 
                      open(runner.args.ccsFofn).read().splitlines()}
        
        # fill in the CCS spot.
        for k,v in zmwsForBCs.items():
            l = []
            for zmw,lZmw in v:
                r = ccsReaders[movieNameFromFile(zmw.baxH5.file.filename)]
                l.append((zmw,lZmw,r[zmw.holeNumber]))
            zmwsForBCs[k] = l
    else:
        # add none to the CCS spot.
        zmwsForBCs = {k:[(zmw,lZmw,None) for zmw,lZmw in v] 
                      for k,v in zmwsForBCs.iteritems()}

    readAndReads = { k:makeReadAndReads(v) for k,v in zmwsForBCs.items() }

    # remove barcodes that don't have a seed read and a set of useable reads.
    readAndReads = { k:v for k,v in readAndReads.items() if v[0] and v[1] }
   
    ## generate FASTA files
    outDir = runner.args.outDir

    for barcode, reads in readAndReads.items():
        bcdir = '/'.join((outDir, barcode))
        if not os.path.exists(bcdir):
            os.makedirs(bcdir)

        ## emit the seeds to separte files
        with FastaWriter("%s/seed_read.fasta" % bcdir) as w:
            w.writeRecord(FastaRecord(reads[0].readName, reads[0].basecalls()))
        
        ## emit the subreads to a single file
        with FastaWriter("%s/subreads.fasta" % bcdir) as w:
            for r in reads[1]:
                w.writeRecord(FastaRecord(r.readName, r.basecalls()))
                
                subreads = reads[1]
        
        ## now, make region files for an eventual blasr alignment
        regs = [(a, h5.File(a, 'r')['/PulseData/Regions']) for a in 
                open(runner.args.inputFofn, 'r').read().splitlines()]
        nfofn = []
        for inFof,regTbl in regs:
            holes = n.in1d(regTbl[:, 0], n.array([ a.holeNumber for a in subreads ]))
            if any(holes): 
                reg = regTbl[holes, :]
            else:
                reg = n.empty(shape = (0, regTbl.shape[1]), dtype = 'int32')
            fname = "%s/%s.rgn.h5" % (bcdir, os.path.basename(inFof))
            nfile = h5.File(fname, 'w')
            ndset = nfile.create_dataset('/PulseData/Regions', data = reg, 
                                         maxshape = (None, None))
            copyAttributes(regTbl, ndset)
            nfile.close()
            nfofn.append(fname)
        
        ofile = open('%s/region.fofn' % bcdir, 'w')
        ofile.writelines("\n".join(nfofn))
        ofile.close()
    
    ## call gcon
    outDirs  = [ (outDir, k) for k in readAndReads.keys() ]
    if runner.args.nProcs == 1:
	outFasta = filter(lambda z: z, map(gconFunc, outDirs))
    else:
    	pool = Pool(runner.args.nProcs)
    	outFasta = filter(lambda z : z, pool.map(gconFunc, outDirs))

    ## write the results
    with FastaWriter('/'.join((outDir, "consensus.fa"))) as w:
        for r in outFasta:
            w.writeRecord(r)

    ## optionally cleanup
    if not runner.args.keepTmpDir:
        for barcode, reads in readAndReads.items():
             bcdir = '/'.join((outDir, barcode))
             shutil.rmtree(bcdir)
        

class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))
        subparsers = self.subParsers
                
        desc = ['Creates a barcode.h5 file from base h5 files.']
        parser_m = subparsers.add_parser('labelZmws', description = "\n".join(desc), 
                                         help = 'Label zmws with barcode annotation',
                                         formatter_class = \
                                             argparse.ArgumentDefaultsHelpFormatter)
        parser_m.add_argument('--outDir', 
                              help = 'Where to write the newly created barcode.h5 files.',
                              default = os.getcwd())
        parser_m.add_argument('--outFofn', help = 'Write to outFofn',
                              default = 'barcode.fofn')
        parser_m.add_argument('--adapterSidePad', help = 'Pad with adapterSidePad bases',
                              default = 4, type = int)
        parser_m.add_argument('--insertSidePad', help = 'Pad with insertSidePad bases',
                              default = 4, type = int)
        parser_m.add_argument('--scoreMode', 
                              help = 'The mode in which the barcodes should be scored.',
                              choices = SCORE_MODES, default = 'symmetric', type = str)
        parser_m.add_argument('--maxAdapters', type = int, default = 20, 
                              help = 'Only score the first maxAdapters')
        parser_m.add_argument('--scoreFirst', action = 'store_true', default = False,
                              help = 'Whether to try to score the leftmost barcode in a trace.')
        parser_m.add_argument('--startTimeCutoff', 
                              help = 'Reads must start before this value in order to be ' + \
                                  'included when scoreFirst is set.', type = float, 
                              default = 10.0)
        parser_m.add_argument('--nZmws', type = int, default = -1, 
                              help = 'Use the first n ZMWs for testing')
        parser_m.add_argument('--nProcs', type = int, default = 8, 
                              help = 'How many processes to use')
        parser_m.add_argument('--saveExtendedInfo', action = 'store_true', default = False,\
                                  help = 'Whether to save extended information to' + \
                                  'the barcode.h5 files; this information is useful for ' + \
                                  'debugging and chimera detection')
        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta', 
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.fofn',
                              help = 'Input base fofn')

        def addFilteringOpts(parser, justBarcode = False):
            ## These are independent of the barcode scoring
            if not justBarcode: 
                parser.add_argument('--minMaxInsertLength', default = 0, type = int, 
                                    help = "ZMW Filter: exclude ZMW if the longest subread" + \
                                        "is less than this amount")
                parser.add_argument('--hqStartTime', default = float("inf"), type = float,
                                    help = "ZMW Filter: exclude ZMW if start time of HQ region" + \
                                        "greater than this value (seconds)")
                parser.add_argument('--minReadScore', default = 0, type = float,
                                    help = "ZMW Filter: exclude ZMW if readScore is less than" + \
                                        "this value")
     
            ## These obviously need the barcode score
            parser.add_argument('--minAvgBarcodeScore', default = 0.0, type = float,
                                help = "ZMW Filter: exclude ZMW if average barcode score " + \
                                    "is less than this value")
            parser.add_argument('--minNumBarcodes', default = 1, type = int,
                                help = "ZMW Filter: exclude ZMW if number of barcodes observed " + \
                                "is less than this value")
            parser.add_argument('--minScoreRatio', default = 1.0, type = float,
                                help = "ZMW Filter: exclude ZMWs whose best score divided by " + \
                                    "the 2nd best score is less than this ratio")

            # Not yet implemented
            # parser.add_argument('--filterChimeras', default = False, action = 'store_true',
            #                     help = "ZMW Filter: exclude ZMWs that appear to be chimeric")

        
        desc = ['Adds information about barcode alignments to a cmp.h5 file',
                'from a previous call to "labelZmws".']
        parser_s = subparsers.add_parser('labelAlignments', description = "\n".join(desc),
                                         help = "Label reads from a barcode or region h5 file",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        addFilteringOpts(parser_s, justBarcode = True)
        parser_s.add_argument('inputFofn', metavar = 'barcode.fofn',
                              help = 'input barcode fofn file')
        parser_s.add_argument('cmpH5', metavar = 'aligned_reads.cmp.h5',
                              help = 'cmp.h5 file to add barcode labels')
       
        desc = ['Takes a bas.h5 fofn and a barcode.h5 fofn and produces',
                'a fast[a|q] file for each barcode.']
        parser_s = subparsers.add_parser('emitFastqs', description = "\n".join(desc),
                                         help = "Write fastq files", 
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser_s.add_argument('--outDir', metavar = 'output.dir',
                              help = 'output directory to write fastq files',
                              default = os.getcwd())

        parser_s.add_argument('--subreads', 
                              help = 'whether to produce fastq files for the subreads;' + \
                                  'the default is to use the CCS reads. This option only' + \
                                  'applies when input.fofn has both consensus and raw reads,' + \
                                  'otherwise the read type from input.fofn will be returned.',
                              action = 'store_true',
                              default = False)
        parser_s.add_argument('--unlabeledZmws', 
                              help = 'whether to emit a fastq file for the unlabeled ZMWs.' + \
                              ' These are the ZMWs where no adapters are found typically',
                              action = 'store_true',
                              default = False)

        parser_s.add_argument('--trim', help = 'trim off barcodes and any excess constant sequence',
                              default = 20, type = int)
        parser_s.add_argument('--fasta', help = ('whether the files produced should be FASTA files as' +
                                                 'opposed to FASTQ'),
                              action = 'store_true',
                              default = False)
        addFilteringOpts(parser_s)
        parser_s.add_argument('inputFofn', metavar = 'input.fofn',
                              help = 'input base or CCS fofn file')
        parser_s.add_argument('barcodeFofn', metavar = 'barcode.fofn',
                              help = 'input barcode.h5 fofn file')

        desc = ['Compute consensus sequences for each barcode.']
        parser_s = subparsers.add_parser('consensus', description = "\n".join(desc),
                                         help = "Compute a consensus sequence for each barcode." + \
                                             "This command relies on the presence of pbdagcon",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser_s.add_argument('--subsample', default = 1, type = float,
                              help = "Subsample ZMWs")
        parser_s.add_argument('--nZmws', default = -1, type = int,
                              help = "Take n ZMWs")
        parser_s.add_argument('--outDir', default = '.', type = str,
                              help = "Use this directory to output results")
        parser_s.add_argument('--keepTmpDir', action = 'store_true', default = False)
        parser_s.add_argument('--ccsFofn', default = '', type = str,
                              help = 'Obtain CCS data from ccsFofn instead of input.fofn')
        parser_s.add_argument('--nProcs', default = 16, type = int,
                              help = 'Use nProcs to execute.')
        addFilteringOpts(parser_s)
        parser_s.add_argument('inputFofn', metavar = 'input.fofn',
                              help = 'input bas.h5 fofn file')
        parser_s.add_argument('barcodeFofn', metavar = 'barcode.fofn',
                              help = 'input bc.h5 fofn file')

    def getVersion(self):
        return __version__

    def run(self):
        logging.debug("Arguments" + str(self.args))
        
        if self.args.subCommand == 'labelZmws':
            makeBarcodeFofnFromBasFofn()
        elif self.args.subCommand == 'labelAlignments':
            labelAlignments()
        elif self.args.subCommand == 'emitFastqs':
            emitFastqs()
        elif self.args.subCommand == 'consensus':
            callConsensus()
        else:
            sys.exit(1)
          
if __name__ == '__main__':    
    runner = Pbbarcode()
    sys.exit(runner.start())
