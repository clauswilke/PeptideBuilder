from PeptideBuilder import Geometry

# test all geometries for correct parameters
def test_geometry_A():
    g = Geometry.geometry("A")
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_C_O_diangle == -60.5
    assert g.N_CA_C_angle == 111.068
    assert g.N_C_CA_CB_diangle == 122.686
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "A"


