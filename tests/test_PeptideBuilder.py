import PeptideBuilder
from Bio.PDB.Polypeptide import PPBuilder


# test whether the correct sequence can be recovered from the structure
def test_sequence_recovery():
    # Build a peptide containing all 20 amino acids
    structure = PeptideBuilder.initialize_res('A')

    for aa in "CDEFGHIKLMNPQRSTVWY":
        structure = PeptideBuilder.add_residue(structure, aa)

    # Extract peptide from structure and compare to expected
    ppb = PPBuilder()
    pp = next(iter(ppb.build_peptides(structure)))
    assert pp.get_sequence() == "ACDEFGHIKLMNPQRSTVWY"
