#!/usr/bin/python

'''
This script builds a variety of structures using all the major functionality in PeptideBuilder. Compare the generated structures to the provided reference structures to see if everything is correct.
'''

from __future__ import print_function
from PeptideBuilder import Geometry
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
outfile = open("test3.txt", 'w')
for aa in "ACDEFGHIKLMNPQRSTVWY":
    print(outfile, Geometry.geometry(aa))

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


# Build a helix containing all 20 amino acids, with slowly varying
# backbone angles, using make_structure().
# The resulting structure should be identical to test2.pdb
phi_list = []
psi_im1_list = []

for i in range(1,20):
    phi_list.append(-60+i)
    psi_im1_list.append(-40-i)
structure = PeptideBuilder.make_structure("ACDEFGHIKLMNPQRSTVWY",\
                                                phi_list, psi_im1_list)
out.set_structure(structure)
out.save("test6.pdb")

# Build a helix containing all 20 amino acids, with slowly varying
# backbone angles, using make_structure().
# The first half of the resulting structure should be identical to
# test6.pdb, while the second half should be slightly different.
phi_list = []
psi_im1_list = []
omega_list = []

for i in range(1,20):
    phi_list.append(-60+i)
    psi_im1_list.append(-40-i)
    omega_list.append(180)
    
for i in range(9,19):
    omega_list[i] = -178
    
structure = PeptideBuilder.make_structure("ACDEFGHIKLMNPQRSTVWY",\
                                                phi_list, psi_im1_list,\
                                                omega_list)
out.set_structure(structure)
out.save("test7.pdb")

