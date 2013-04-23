#!/usr/bin/python

'''
This script builds a variety of structures using all the major functionality in PeptideBuilder. Compare the generated structures to the provided reference structures to see if everything is correct.
'''


import Geometry
import PeptideBuilder
import Bio.PDB


# Build a peptide containing all 20 amino acids
structure = PeptideBuilder.initialize_res('A')

for aa in "CDEFGHIKLMNPQRSTVWY":
    structure = PeptideBuilder.add_residue(structure, aa)
    
out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("test1.pdb")

# Build a helix containing all 20 amino acids, with slowly varying backbone angles
phi = -60
psi_im1 = -40
geo = Geometry.geometry('A')
geo.phi = phi
geo.psi_im1 = psi_im1
structure = PeptideBuilder.initialize_res(geo)

for aa in "CDEFGHIKLMNPQRSTVWY":
    phi += 1
    psi_im1 -= 1
    geo = Geometry.geometry(aa)
    geo.phi = phi
    geo.psi_im1 = psi_im1
    structure = PeptideBuilder.add_residue(structure, geo)
    
out.set_structure(structure)
out.save("test2.pdb")

# Print out all geometries
outfile = file("test3.txt", 'w')
for aa in "ACDEFGHIKLMNPQRSTVWY":
    print >>outfile, Geometry.geometry(aa)

# Build a helix containing all 20 amino acids from list of geometries.
# The structure should be identical to test1.pdb

geos = []
for aa in "ACDEFGHIKLMNPQRSTVWY":
    geos.append(Geometry.geometry(aa))
structure = PeptideBuilder.make_structure_from_geos(geos)
out.set_structure(structure)
out.save("test4.pdb")

# Build a peptide containing all 20 amino acids in extended conformation.
# The structure should be identical to test1.pdb

structure = PeptideBuilder.make_extended_structure("ACDEFGHIKLMNPQRSTVWY")
out.set_structure(structure)
out.save("test4.pdb")

# Build a peptide containing all 20 amino acids from list of geometries.
# The structure should be identical to test1.pdb

geos = []
for aa in "ACDEFGHIKLMNPQRSTVWY":
    geos.append(Geometry.geometry(aa))
structure = PeptideBuilder.make_structure_from_geos(geos)  
out.set_structure(structure)
out.save("test5.pdb")