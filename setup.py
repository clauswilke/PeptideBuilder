#!/usr/bin/python

##############################################################################
##  PeptideBuilder: builds peptide models in PDB format using geometrical inputs
##
##  Written by Matthew Z Tien (Matthew.Tien89@gmail.com) 
##############################################################################

'''
Setup.py script (uses setuptools) for building, testing, and installing PeptideBuilder.
To build and install the package as root (globally), enter (from this directory!) - 
    sudo python setup.py build
    sudo python setup.py install
    
To install for a particular user (locally), enter - 
    python setup.py build
    python setup.py build --user # where user is the computer account to install PeptideBuilder in
'''



from setuptools import setup
setup(name = 'PeptideBuilder', 
    version = '1.3', 
    description = 'Tools to create peptide PDB files using geometry as input',
    author = 'Matthew Z. Tien', 
    author_email = 'Matthew.Tien89@gmail.com', 
    url = 'https://github.com/mtien/PeptideBuilder',
    download_url = 'https://github.com/smtien/PeptideBuilder/mypackage/tarball/1.3',
    platforms = 'Tested on Mac OS X.',
    packages = ['PeptideBuilder'],
    install_requires=['Bio.PDB']
)
