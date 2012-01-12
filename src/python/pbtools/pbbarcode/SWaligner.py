from ctypes import CDLL
import os
import numpy
import pkg_resources

class SWaligner(object):
    def __init__(self):
        # setup.py should put sw.so in the following path.
        self.SW_DLL_PATH = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "sw.so" 
        self._dll        = CDLL(self.SW_DLL_PATH)
        self.dpMat       = self._dll.allocate_dp_mat()

    def score(self, tSeq, qSeq):
        return self._dll.compute_align_score(self.dpMat, tSeq, qSeq)
