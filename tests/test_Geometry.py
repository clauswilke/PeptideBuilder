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


def test_geometry_C():
    g = Geometry.geometry("C")
    assert g.CA_CB_SG_angle == 113.8169
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_SG_length == 1.808
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_SG_diangle == -62.2
    assert g.N_CA_C_O_diangle == -60.0
    assert g.N_CA_C_angle == 110.8856
    assert g.N_C_CA_CB_diangle == 122.5037
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "C"


def test_geometry_D():
    g = Geometry.geometry("D")
    assert g.CA_CB_CG_OD1_diangle == -46.7
    assert g.CA_CB_CG_OD2_diangle == 133.3
    assert g.CA_CB_CG_angle == 113.06
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.51
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_OD1_angle == 119.22
    assert g.CB_CG_OD2_angle == 118.218
    assert g.CB_CG_length == 1.52
    assert g.CG_OD1_length == 1.25
    assert g.CG_OD2_length == 1.25
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -66.4
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 111.03
    assert g.N_C_CA_CB_diangle == 122.82
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "D"


def test_geometry_E():
    g = Geometry.geometry("E")
    assert g.CA_CB_CG_CD_diangle == -179.8
    assert g.CA_CB_CG_angle == 113.82
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.511
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD_OE1_diangle == -6.2
    assert g.CB_CG_CD_OE2_diangle == 173.8
    assert g.CB_CG_CD_angle == 113.31
    assert g.CB_CG_length == 1.52
    assert g.CD_OE1_length == 1.25
    assert g.CD_OE2_length == 1.25
    assert g.CG_CD_OE1_angle == 119.02
    assert g.CG_CD_OE2_angle == 118.08
    assert g.CG_CD_length == 1.52
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -63.8
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 111.1703
    assert g.N_C_CA_CB_diangle == 122.8702
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "E"


def test_geometry_F():
    g = Geometry.geometry("F")
    assert g.CA_CB_CG_CD1_diangle == 93.3
    assert g.CA_CB_CG_CD2_diangle == -86.7
    assert g.CA_CB_CG_angle == 113.85
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5316
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD1_CE1_diangle == 180.0
    assert g.CB_CG_CD1_angle == 120.0
    assert g.CB_CG_CD2_CE2_diangle == 180.0
    assert g.CB_CG_CD2_angle == 120.0
    assert g.CB_CG_length == 1.5
    assert g.CD1_CE1_CZ_angle == 120.0
    assert g.CD1_CE1_length == 1.39
    assert g.CD2_CE2_length == 1.39
    assert g.CE1_CZ_length == 1.39
    assert g.CG_CD1_CE1_CZ_diangle == 0.0
    assert g.CG_CD1_CE1_angle == 120.0
    assert g.CG_CD1_length == 1.39
    assert g.CG_CD2_CE2_angle == 120.0
    assert g.CG_CD2_length == 1.39
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -64.7
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.7528
    assert g.N_C_CA_CB_diangle == 122.6054
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "F"


def test_geometry_G():
    g = Geometry.geometry("G")
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5117
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_C_O_diangle == 180.0
    assert g.N_CA_C_angle == 110.8914
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "G"


def test_geometry_H():
    g = Geometry.geometry("H")
    assert g.CA_CB_CG_CD2_diangle == 104.3
    assert g.CA_CB_CG_ND1_diangle == -75.7
    assert g.CA_CB_CG_angle == 113.74
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.4732
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD2_NE2_diangle == 180.0
    assert g.CB_CG_CD2_angle == 130.61
    assert g.CB_CG_ND1_CE1_diangle == 180.0
    assert g.CB_CG_ND1_angle == 122.85
    assert g.CB_CG_length == 1.49
    assert g.CD2_NE2_length == 1.35
    assert g.CG_CD2_NE2_angle == 108.5
    assert g.CG_CD2_length == 1.35
    assert g.CG_ND1_CE1_angle == 108.5
    assert g.CG_ND1_length == 1.38
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.ND1_CE1_length == 1.32
    assert g.N_CA_CB_CG_diangle == -63.2
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 111.0859
    assert g.N_C_CA_CB_diangle == 122.6711
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "H"


def test_geometry_I():
    g = Geometry.geometry("I")
    assert g.CA_CB_CG1_CD1_diangle == 169.8
    assert g.CA_CB_CG1_angle == 110.7
    assert g.CA_CB_CG2_angle == 110.4
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5403
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG1_CD1_angle == 113.97
    assert g.CB_CG1_length == 1.527
    assert g.CB_CG2_length == 1.527
    assert g.CG1_CD1_length == 1.52
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG1_diangle == 59.7
    assert g.N_CA_CB_CG2_diangle == -61.6
    assert g.N_CA_C_O_diangle == -60.0
    assert g.N_CA_C_angle == 109.7202
    assert g.N_C_CA_CB_diangle == 123.2347
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "I"


def test_geometry_K():
    g = Geometry.geometry("K")
    assert g.CA_CB_CG_CD_diangle == -178.1
    assert g.CA_CB_CG_angle == 113.83
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.54
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD_CE_diangle == -179.6
    assert g.CB_CG_CD_angle == 111.79
    assert g.CB_CG_length == 1.52
    assert g.CD_CE_NZ_angle == 124.79
    assert g.CD_CE_length == 1.46
    assert g.CE_NZ_length == 1.33
    assert g.CG_CD_CE_NZ_diangle == 179.6
    assert g.CG_CD_CE_angle == 111.68
    assert g.CG_CD_length == 1.52
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -64.5
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 111.08
    assert g.N_C_CA_CB_diangle == 122.76
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "K"


def test_geometry_L():
    g = Geometry.geometry("L")
    assert g.CA_CB_CG_CD1_diangle == 174.9
    assert g.CA_CB_CG_CD2_diangle == 66.7
    assert g.CA_CB_CG_angle == 116.1
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.4647
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD1_angle == 110.27
    assert g.CB_CG_CD2_angle == 110.58
    assert g.CB_CG_length == 1.53
    assert g.CG_CD1_length == 1.524
    assert g.CG_CD2_length == 1.525
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -60.1
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.8652
    assert g.N_C_CA_CB_diangle == 122.4948
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "L"


def test_geometry_M():
    g = Geometry.geometry("M")
    assert g.CA_CB_CG_SD_diangle == -179.6
    assert g.CA_CB_CG_angle == 113.68
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.4816
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_SD_CE_diangle == 70.1
    assert g.CB_CG_SD_angle == 112.69
    assert g.CB_CG_length == 1.52
    assert g.CG_SD_CE_angle == 100.61
    assert g.CG_SD_length == 1.81
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -64.4
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.9416
    assert g.N_C_CA_CB_diangle == 122.6733
    assert g.SD_CE_length == 1.79
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "M"


def test_geometry_N():
    g = Geometry.geometry("N")
    assert g.CA_CB_CG_ND2_diangle == 121.7
    assert g.CA_CB_CG_OD1_diangle == -58.3
    assert g.CA_CB_CG_angle == 112.62
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.4826
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_ND2_angle == 116.48
    assert g.CB_CG_OD1_angle == 120.85
    assert g.CB_CG_length == 1.52
    assert g.CG_ND2_length == 1.33
    assert g.CG_OD1_length == 1.23
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -65.5
    assert g.N_CA_C_O_diangle == -60.0
    assert g.N_CA_C_angle == 111.5
    assert g.N_C_CA_CB_diangle == 123.2254
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "N"


def test_geometry_P():
    g = Geometry.geometry("P")
    assert g.CA_CB_CG_CD_diangle == -34.8
    assert g.CA_CB_CG_angle == 104.21
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.2945
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD_angle == 105.03
    assert g.CB_CG_length == 1.49
    assert g.CG_CD_length == 1.5
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == 29.6
    assert g.N_CA_C_O_diangle == -45.0
    assert g.N_CA_C_angle == 112.7499
    assert g.N_C_CA_CB_diangle == 115.2975
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "P"


def test_geometry_Q():
    g = Geometry.geometry("Q")
    assert g.CA_CB_CG_CD_diangle == -69.6
    assert g.CA_CB_CG_angle == 113.75
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5029
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD_NE2_diangle == 129.5
    assert g.CB_CG_CD_OE1_diangle == -50.5
    assert g.CB_CG_CD_angle == 112.78
    assert g.CB_CG_length == 1.52
    assert g.CD_NE2_length == 1.33
    assert g.CD_OE1_length == 1.24
    assert g.CG_CD_NE2_angle == 116.5
    assert g.CG_CD_OE1_angle == 120.86
    assert g.CG_CD_length == 1.52
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -60.2
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 111.0849
    assert g.N_C_CA_CB_diangle == 122.8134
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "Q"


def test_geometry_R():
    g = Geometry.geometry("R")
    assert g.CA_CB_CG_CD_diangle == -179.2
    assert g.CA_CB_CG_angle == 113.83
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.54
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD_NE_diangle == -179.3
    assert g.CB_CG_CD_angle == 111.79
    assert g.CB_CG_length == 1.52
    assert g.CD_NE_CZ_NH1_diangle == 0.0
    assert g.CD_NE_CZ_NH2_diangle == 180.0
    assert g.CD_NE_CZ_angle == 124.79
    assert g.CD_NE_length == 1.46
    assert g.CG_CD_NE_CZ_diangle == -178.7
    assert g.CG_CD_NE_angle == 111.68
    assert g.CG_CD_length == 1.52
    assert g.CZ_NH1_length == 1.33
    assert g.CZ_NH2_length == 1.33
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.NE_CZ_NH1_angle == 120.64
    assert g.NE_CZ_NH2_angle == 119.63
    assert g.NE_CZ_length == 1.33
    assert g.N_CA_CB_CG_diangle == -65.2
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.98
    assert g.N_C_CA_CB_diangle == 122.76
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "R"


def test_geometry_S():
    g = Geometry.geometry("S")
    assert g.CA_CB_OG_angle == 110.773
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_OG_length == 1.417
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_OG_diangle == -63.3
    assert g.N_CA_C_O_diangle == -60.0
    assert g.N_CA_C_angle == 111.2812
    assert g.N_C_CA_CB_diangle == 122.6618
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "S"


def test_geometry_T():
    g = Geometry.geometry("T")
    assert g.CA_CB_CG2_angle == 111.13
    assert g.CA_CB_OG1_angle == 109.18
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5359
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG2_length == 1.53
    assert g.CB_OG1_length == 1.43
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG2_diangle == -60.3
    assert g.N_CA_CB_OG1_diangle == 60.0
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.7014
    assert g.N_C_CA_CB_diangle == 123.0953
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "T"


def test_geometry_V():
    g = Geometry.geometry("V")
    assert g.CA_CB_CG1_angle == 110.7
    assert g.CA_CB_CG2_angle == 110.4
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5686
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG1_length == 1.527
    assert g.CB_CG2_length == 1.527
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG1_diangle == 177.2
    assert g.N_CA_CB_CG2_diangle == -63.3
    assert g.N_CA_C_O_diangle == -60.0
    assert g.N_CA_C_angle == 109.7698
    assert g.N_C_CA_CB_diangle == 123.2347
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "V"


def test_geometry_W():
    g = Geometry.geometry("W")
    assert g.CA_CB_CG_CD1_diangle == 96.3
    assert g.CA_CB_CG_CD2_diangle == -83.7
    assert g.CA_CB_CG_angle == 114.1
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5117
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD1_NE1_diangle == 180.0
    assert g.CB_CG_CD1_angle == 127.07
    assert g.CB_CG_CD2_CE2_diangle == 180.0
    assert g.CB_CG_CD2_CE3_diangle == 0.0
    assert g.CB_CG_CD2_angle == 126.66
    assert g.CB_CG_length == 1.5
    assert g.CD1_NE1_length == 1.38
    assert g.CD2_CE2_CZ2_CH2_diangle == 0.0
    assert g.CD2_CE2_CZ2_angle == 120.0
    assert g.CD2_CE2_length == 1.4
    assert g.CD2_CE3_CZ3_angle == 120.0
    assert g.CD2_CE3_length == 1.4
    assert g.CE2_CZ2_CH2_angle == 120.0
    assert g.CE2_CZ2_length == 1.4
    assert g.CE3_CZ3_length == 1.4
    assert g.CG_CD1_NE1_angle == 108.5
    assert g.CG_CD1_length == 1.37
    assert g.CG_CD2_CE2_CZ2_diangle == 180.0
    assert g.CG_CD2_CE2_angle == 108.5
    assert g.CG_CD2_CE3_CZ3_diangle == 180.0
    assert g.CG_CD2_CE3_angle == 133.83
    assert g.CG_CD2_length == 1.43
    assert g.CZ2_CH2_length == 1.4
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -66.4
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.8914
    assert g.N_C_CA_CB_diangle == 122.6112
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "W"


def test_geometry_Y():
    g = Geometry.geometry("Y")
    assert g.CA_CB_CG_CD1_diangle == 93.1
    assert g.CA_CB_CG_CD2_diangle == 273.1
    assert g.CA_CB_CG_angle == 113.8
    assert g.CA_CB_length == 1.52
    assert g.CA_C_N_angle == 116.642992978143
    assert g.CA_C_O_angle == 120.5434
    assert g.CA_C_length == 1.52
    assert g.CA_N_length == 1.46
    assert g.CB_CG_CD1_CE1_diangle == 180.0
    assert g.CB_CG_CD1_angle == 120.98
    assert g.CB_CG_CD2_CE2_diangle == 180.0
    assert g.CB_CG_CD2_angle == 120.82
    assert g.CB_CG_length == 1.51
    assert g.CD1_CE1_CZ_OH_diangle == 180.0
    assert g.CD1_CE1_CZ_angle == 120.0
    assert g.CD1_CE1_length == 1.39
    assert g.CD2_CE2_length == 1.39
    assert g.CE1_CZ_OH_angle == 119.78
    assert g.CE1_CZ_length == 1.39
    assert g.CG_CD1_CE1_CZ_diangle == 0.0
    assert g.CG_CD1_CE1_angle == 120.0
    assert g.CG_CD1_length == 1.39
    assert g.CG_CD2_CE2_angle == 120.0
    assert g.CG_CD2_length == 1.39
    assert g.CZ_OH_length == 1.39
    assert g.C_CA_CB_angle == 109.5
    assert g.C_N_CA_angle == 121.382215820277
    assert g.C_O_length == 1.23
    assert g.N_CA_CB_CG_diangle == -64.3
    assert g.N_CA_C_O_diangle == 120.0
    assert g.N_CA_C_angle == 110.9288
    assert g.N_C_CA_CB_diangle == 122.6023
    assert g.omega == 180.0
    assert g.peptide_bond == 1.33
    assert g.phi == -120
    assert g.psi_im1 == 140
    assert g.residue_name == "Y"
