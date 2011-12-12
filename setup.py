from setuptools import setup, Extension, find_packages

setup(
    name = 'pbtools.barcode',
    version='0.1.0',
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    scripts = ['src/barcodeLabeler.py'],
    packages = find_packages('src'),  
    package_dir = {'':'src'},
    namespace_packages = ['pbtools'],
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.1',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0'
        ]
    )
