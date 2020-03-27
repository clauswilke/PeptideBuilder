import PeptideBuilder
from PeptideBuilder import Geometry
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBParser
import os

def compare_residues(r1, r2):
    result = True
    result = result and r1 == r2
    for a1, a2 in zip(r1, r2):
        #print(a1.coord - a2.coord)
        result = result and (abs(a1.coord - a2.coord) < 0.001).all()
    return result

def compare_to_reference(structure, ref_file):
    file = os.path.join("tests", "pdbs", ref_file)
    parser = PDBParser()
    ref_structure = parser.get_structure('test', file)
    
    result = True
    res = list(list(structure[0])[0])
    ref_res = list(list(ref_structure[0])[0])
    for r1, r2 in zip(res, ref_res):
        result = result and compare_residues(r1, r2)
    return result

# Build a peptide containing all 20 amino acids
def test_add_residue():
    structure = PeptideBuilder.initialize_res('A')
    for aa in "CDEFGHIKLMNPQRSTVWY":
        structure = PeptideBuilder.add_residue(structure, aa)

    # extract peptide from structure and compare to expected
    ppb = PPBuilder()
    pp = next(iter(ppb.build_peptides(structure)))
    assert pp.get_sequence() == "ACDEFGHIKLMNPQRSTVWY"

    # now compare to saved reference structure
    assert compare_to_reference(structure, "extended.pdb")
        
# Build a helix containing all 20 amino acids, with slowly varying backbone angles
def test_add_residue2():
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
        
    # now compare to saved reference structure
    assert compare_to_reference(structure, "helix.pdb")


# Build a helix containing all 20 amino acids from list of geometries.
# The structure should be identical to `extended.pdb`
def test_make_structure_from_geos():
    geos = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        geos.append(Geometry.geometry(aa))
    structure = PeptideBuilder.make_structure_from_geos(geos)
    assert compare_to_reference(structure, "extended.pdb")



# Build a peptide containing all 20 amino acids in extended conformation.
# The structure should be identical to `extended.pdb`
def test_make_extended_structure():
    structure = PeptideBuilder.make_extended_structure("ACDEFGHIKLMNPQRSTVWY")
    assert compare_to_reference(structure, "extended.pdb")


# Build a peptide containing all 20 amino acids from list of geometries.
# The structure should be identical to `extended.pdb`
def test_make_structure_from_geos2():
    geos = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        geos.append(Geometry.geometry(aa))
    structure = PeptideBuilder.make_structure_from_geos(geos)
    assert compare_to_reference(structure, "extended.pdb")


# Build a helix containing all 20 amino acids, with slowly varying
# backbone angles, using make_structure().
# The resulting structure should be identical to `helix.pdb`
def test_make_structure():
    phi_list = []
    psi_im1_list = []

    for i in range(1, 20):
        phi_list.append(-60+i)
        psi_im1_list.append(-40-i)
    structure = PeptideBuilder.make_structure(
        "ACDEFGHIKLMNPQRSTVWY",
        phi_list, psi_im1_list
    )
    assert compare_to_reference(structure, "helix.pdb")


# Build a helix containing all 20 amino acids, with slowly varying
# backbone angles, using make_structure(). Now we're changing omega also.
# The first half of the resulting structure should be identical to
# `helix.pdb`, while the second half should be slightly different.
def test_make_structure2():
    phi_list = []
    psi_im1_list = []
    omega_list = []

    for i in range(1,20):
        phi_list.append(-60+i)
        psi_im1_list.append(-40-i)
        omega_list.append(180)
    
    for i in range(9,19):
        omega_list[i] = -178
    
    structure = PeptideBuilder.make_structure(
        "ACDEFGHIKLMNPQRSTVWY",
        phi_list, psi_im1_list, omega_list
    )
    assert compare_to_reference(structure, "helix2.pdb")
