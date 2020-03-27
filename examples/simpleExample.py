#!/usr/bin/python

'''
Simple example script demonstrating how to use the PeptideBuilder library.

The script generates a peptide consisting of six arginines in alpha-helix
conformation, and it stores the peptide under the name "example.pdb".
'''

from __future__ import print_function
from PeptideBuilder import Geometry
import PeptideBuilder


geo = Geometry.geometry('G')
geo.phi=-60
geo.psi_im1=-40
structure = PeptideBuilder.initialize_res(geo)
for i in range(5):
    structure = PeptideBuilder.add_residue(structure, geo)
    
import Bio.PDB
out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save( "example.pdb" )
