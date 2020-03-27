import PeptideBuilder
from PeptideBuilder import Geometry
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBParser

def compare_residues(r1, r2):
    result = True
    result = result and r1 == r2
    for a1, a2 in zip(r1, r2):
        #print(a1.coord - a2.coord)
        result = result and (abs(a1.coord - a2.coord) < 0.001).all()
    return result

def compare_to_reference(structure, ref_file):
    parser = PDBParser()
    ref_structure = parser.get_structure('test', ref_file)
    
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
    assert compare_to_reference(structure, "tests/pdbs/test1.pdb")
        
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
    # now compare to saved reference structure
    assert compare_to_reference(structure, "tests/pdbs/test2.pdb")


# Build a helix containing all 20 amino acids from list of geometries.
# The structure should be identical to test1.pdb
def test_make_structure_from_geos():
    geos = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        geos.append(Geometry.geometry(aa))
    structure = PeptideBuilder.make_structure_from_geos(geos)

    assert compare_to_reference(structure, "tests/pdbs/test1.pdb")



# Build a peptide containing all 20 amino acids in extended conformation.
# The structure should be identical to test1.pdb
def test_make_extended_structure():
    structure = PeptideBuilder.make_extended_structure("ACDEFGHIKLMNPQRSTVWY")

    assert compare_to_reference(structure, "tests/pdbs/test1.pdb")


# Build a peptide containing all 20 amino acids from list of geometries.
# The structure should be identical to test1.pdb
def test_make_structure_from_geos2():
    geos = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        geos.append(Geometry.geometry(aa))
    structure = PeptideBuilder.make_structure_from_geos(geos)  

    assert compare_to_reference(structure, "tests/pdbs/test1.pdb")

