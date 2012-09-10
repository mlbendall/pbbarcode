from setuptools import setup, Extension, find_packages
import os
import sys

#if os.system("make -C src/C all"):
#    sys.exit(1)
    
setup(
    name = 'pbtools.barcode',
    version='0.1.0',
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    scripts = ['src/python/pbbarcode.py'],
    packages = find_packages('src/python'),  
    package_dir = {'':'src/python'},
    namespace_packages = ['pbtools'],
    ext_modules=[Extension('pbtools/pbbarcode/sw', ['src/C/sw.c'], extra_compile_args=["-O3","-shared"])], 
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.1',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0'
        ]
    )
