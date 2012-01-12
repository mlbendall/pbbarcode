from ctypes import *
import os
import numpy
import pkg_resources

# setup.py should put sw.so in the following path.
SW_DLL_PATH = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "sw.so" 

