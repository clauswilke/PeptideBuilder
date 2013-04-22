#!/usr/bin/python

'''
This script builds a variety of structures using all the major functionality in PeptideBuilder. Compare the generated structures to the provided reference structures to see if everything is correct.
'''


import Geometry
import PeptideBuilder
import Bio.PDB


# Build a helix containing all 20 amino acids
structure = PeptideBuilder.initialize_res('A')

for aa in "CDEFGHIKLMNPQRSTVWY":
	structure = PeptideBuilder.add_residue(structure, aa)
    
out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save( "test1.pdb" )

# Build a helix containing all 20 amino acids, with slowly varying backbone angles
phi = -60
psi_im1 = -40
geo = Geometry.geometry('A')
#geo.phi = phi
#geo.psi_im1 = psi_im1
structure = PeptideBuilder.initialize_res(geo)

for aa in "CDEFGHIKLMNPQRSTVWY":
	phi += 1
	psi_im1 -= 1
	geo = Geometry.geometry(aa)
#	geo.phi = phi
#	geo.psi_im1 = psi_im1
	structure = PeptideBuilder.add_residue(structure, geo)
    
out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save( "test2.pdb" )


