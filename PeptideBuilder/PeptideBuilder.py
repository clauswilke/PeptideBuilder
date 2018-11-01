'''This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The PeptideBuilder module contains code to generate 3D
structures of peptides. It requires the Geometry module
(also part of the PeptideBuilder library), which contains
default bond lengths and angles for all amino acids.

This module also requires the Bio.PDB module from
Biopython, for structure manipulation.

This file is provided to you under the MIT License.'''

from __future__ import print_function
from Bio.PDB import *
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Structure import *
from Bio.PDB.Vector import *
from Bio.PDB.Entity import*
from .Geometry import *
import math, warnings
import numpy


def get_prop(atm):
    print(atm.get_name())
    print(atm.get_coord())
    print(atm.get_vector())
    print(atm.get_bfactor())
    print(atm.get_anisou())
    print(atm.get_occupancy())
    print(atm.get_altloc())
    print(atm.get_fullname())
    print(atm.get_serial_number())
    print(atm.get_parent())
    print(atm.get_id())
    print(atm.get_full_id())
    print(atm.get_level())

def calculateCoordinates(refA, refB, refC, L, ang, di):
    AV=refA.get_vector()
    BV=refB.get_vector()
    CV=refC.get_vector()
    
    CA=AV-CV
    CB=BV-CV

    ##CA vector
    AX=CA[0]
    AY=CA[1]
    AZ=CA[2]

    ##CB vector
    BX=CB[0]
    BY=CB[1]
    BZ=CB[2]

    ##Plane Parameters
    A=(AY*BZ)-(AZ*BY)
    B=(AZ*BX)-(AX*BZ)
    G=(AX*BY)-(AY*BX)

    ##Dot Product Constant
    F= math.sqrt(BX*BX + BY*BY + BZ*BZ) * L * math.cos(ang*(math.pi/180.0))

    ##Constants
    const=math.sqrt( math.pow((B*BZ-BY*G),2) *(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L))
    denom= (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G)

    X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+const)/denom

    if((B==0 or BZ==0) and (BY==0 or G==0)):
        const1=math.sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)))
        Y= ((-A*B*X)+const1)/(B*B+G*G)
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G))
    else:
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*math.pow(B*BZ-BY*G,2) + BX*const) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*const)) / ((B*BZ-BY*G)*denom)
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*math.pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*const + A*BY*const) / ((B*BZ-BY*G)*denom)

    
    #GET THE NEW VECTOR from the orgin
    D=Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp=calc_dihedral(AV, BV, CV, D)*(180.0/math.pi)
    
  
    di=di-temp
    rot= rotaxis(math.pi*(di/180.0), CV-BV)
    D=(D-BV).left_multiply(rot)+BV
    
    return D.get_array()

def makeGly(segID, N, CA, C, O, geo):
    '''Creates a Glycine residue'''
    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "GLY", '    ')

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)

    ##print(res)
    return res

def makeAla(segID, N, CA, C, O, geo):
    '''Creates an Alanine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")

    ##Create Residue Data Structure
    res = Residue((' ', segID, ' '), "ALA", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    return res

def makeSer(segID, N, CA, C, O, geo):
    '''Creates a Serine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_OG_length=geo.CB_OG_length
    CA_CB_OG_angle=geo.CA_CB_OG_angle
    N_CA_CB_OG_diangle=geo.N_CA_CB_OG_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    oxygen_g= calculateCoordinates(N, CA, CB, CB_OG_length, CA_CB_OG_angle, N_CA_CB_OG_diangle)
    OG= Atom("OG", oxygen_g, 0.0, 1.0, " ", " OG", 0, "O")

    ##Create Reside Data Structure
    res= Residue((' ', segID, ' '), "SER", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(OG)

    ##print(res)
    return res

def makeCys(segID, N, CA, C, O, geo):
    '''Creates a Cysteine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_SG_length= geo.CB_SG_length
    CA_CB_SG_angle= geo.CA_CB_SG_angle
    N_CA_CB_SG_diangle= geo.N_CA_CB_SG_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    sulfur_g= calculateCoordinates(N, CA, CB, CB_SG_length, CA_CB_SG_angle, N_CA_CB_SG_diangle)
    SG= Atom("SG", sulfur_g, 0.0, 1.0, " ", " SG", 0, "S")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "CYS", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(SG)
    return res

def makeVal(segID, N, CA, C, O, geo):
    '''Creates a Valine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG1_length=geo.CB_CG1_length
    CA_CB_CG1_angle=geo.CA_CB_CG1_angle
    N_CA_CB_CG1_diangle=geo.N_CA_CB_CG1_diangle
    
    CB_CG2_length=geo.CB_CG2_length
    CA_CB_CG2_angle=geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle=geo.N_CA_CB_CG2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g1= calculateCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle)
    CG1= Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle)
    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "VAL", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG1)
    res.add(CG2)
    return res

def makeIle(segID, N, CA, C, O, geo):
    '''Creates an Isoleucine residue'''
    ##R-group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG1_length=geo.CB_CG1_length
    CA_CB_CG1_angle=geo.CA_CB_CG1_angle
    N_CA_CB_CG1_diangle=geo.N_CA_CB_CG1_diangle 
    
    CB_CG2_length=geo.CB_CG2_length
    CA_CB_CG2_angle=geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle= geo.N_CA_CB_CG2_diangle

    CG1_CD1_length= geo.CG1_CD1_length
    CB_CG1_CD1_angle= geo.CB_CG1_CD1_angle
    CA_CB_CG1_CD1_diangle= geo.CA_CB_CG1_CD1_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g1= calculateCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle)
    CG1= Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle)
    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d1= calculateCoordinates(CA, CB, CG1, CG1_CD1_length, CB_CG1_CD1_angle, CA_CB_CG1_CD1_diangle)
    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "ILE", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG1)
    res.add(CG2)
    res.add(CD1)
    return res

def makeLeu(segID, N, CA, C, O, geo):
    '''Creates a Leucine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle= geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_CD1_length=geo.CG_CD1_length
    CB_CG_CD1_angle=geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle=geo.CA_CB_CG_CD1_diangle

    CG_CD2_length=geo.CG_CD2_length
    CB_CG_CD2_angle=geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle=geo.CA_CB_CG_CD2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g1= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g1, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle)
    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle)
    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "LEU", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CD2)
    return res
    
def makeThr(segID, N, CA, C, O, geo):
    '''Creates a Threonine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_OG1_length=geo.CB_OG1_length
    CA_CB_OG1_angle=geo.CA_CB_OG1_angle
    N_CA_CB_OG1_diangle=geo.N_CA_CB_OG1_diangle 
        
    CB_CG2_length=geo.CB_CG2_length
    CA_CB_CG2_angle=geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle= geo.N_CA_CB_CG2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    oxygen_g1= calculateCoordinates(N, CA, CB, CB_OG1_length, CA_CB_OG1_angle, N_CA_CB_OG1_diangle)
    OG1= Atom("OG1", oxygen_g1, 0.0, 1.0, " ", " OG1", 0, "O")
    carbon_g2= calculateCoordinates(N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle)
    CG2= Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "THR", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(OG1)
    res.add(CG2)
    return res

def makeArg(segID, N, CA, C, O, geo):
    '''Creates an Arginie residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle= geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_CD_length=geo.CG_CD_length
    CB_CG_CD_angle=geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle=geo.CA_CB_CG_CD_diangle
    
    CD_NE_length=geo.CD_NE_length
    CG_CD_NE_angle=geo.CG_CD_NE_angle
    CB_CG_CD_NE_diangle=geo.CB_CG_CD_NE_diangle

    NE_CZ_length=geo.NE_CZ_length
    CD_NE_CZ_angle=geo.CD_NE_CZ_angle
    CG_CD_NE_CZ_diangle=geo.CG_CD_NE_CZ_diangle

    CZ_NH1_length=geo.CZ_NH1_length
    NE_CZ_NH1_angle=geo.NE_CZ_NH1_angle
    CD_NE_CZ_NH1_diangle=geo.CD_NE_CZ_NH1_diangle

    CZ_NH2_length=geo.CZ_NH2_length
    NE_CZ_NH2_angle=geo.NE_CZ_NH2_angle
    CD_NE_CZ_NH2_diangle=geo.CD_NE_CZ_NH2_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle)
    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    nitrogen_e= calculateCoordinates(CB, CG, CD, CD_NE_length, CG_CD_NE_angle, CB_CG_CD_NE_diangle)
    NE= Atom("NE", nitrogen_e, 0.0, 1.0, " ", " NE", 0, "N")
    carbon_z= calculateCoordinates(CG, CD, NE, NE_CZ_length, CD_NE_CZ_angle, CG_CD_NE_CZ_diangle)
    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    nitrogen_h1= calculateCoordinates(CD, NE, CZ, CZ_NH1_length, NE_CZ_NH1_angle, CD_NE_CZ_NH1_diangle)
    NH1= Atom("NH1", nitrogen_h1, 0.0, 1.0, " ", " NH1", 0, "N")
    nitrogen_h2= calculateCoordinates(CD, NE, CZ, CZ_NH2_length, NE_CZ_NH2_angle, CD_NE_CZ_NH2_diangle)
    NH2= Atom("NH2", nitrogen_h2, 0.0, 1.0, " ", " NH2", 0, "N")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "ARG", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(NE)
    res.add(CZ)
    res.add(NH1)
    res.add(NH2)
    return res

def makeLys(segID, N, CA, C, O, geo):
    '''Creates a Lysine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD_length=geo.CG_CD_length
    CB_CG_CD_angle=geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle=geo.CA_CB_CG_CD_diangle

    CD_CE_length=geo.CD_CE_length
    CG_CD_CE_angle=geo.CG_CD_CE_angle
    CB_CG_CD_CE_diangle=geo.CB_CG_CD_CE_diangle

    CE_NZ_length=geo.CE_NZ_length
    CD_CE_NZ_angle=geo.CD_CE_NZ_angle
    CG_CD_CE_NZ_diangle=geo.CG_CD_CE_NZ_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle)
    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    carbon_e= calculateCoordinates(CB, CG, CD, CD_CE_length, CG_CD_CE_angle, CB_CG_CD_CE_diangle)
    CE= Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")
    nitrogen_z= calculateCoordinates(CG, CD, CE, CE_NZ_length, CD_CE_NZ_angle, CG_CD_CE_NZ_diangle)
    NZ= Atom("NZ", nitrogen_z, 0.0, 1.0, " ", " NZ", 0, "N")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "LYS", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(CE)
    res.add(NZ)
    return res

def makeAsp(segID, N, CA, C, O, geo):
    '''Creates an Aspartic Acid residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_OD1_length=geo.CG_OD1_length
    CB_CG_OD1_angle=geo.CB_CG_OD1_angle
    CA_CB_CG_OD1_diangle=geo.CA_CB_CG_OD1_diangle

    CG_OD2_length=geo.CG_OD2_length
    CB_CG_OD2_angle=geo.CB_CG_OD2_angle
    CA_CB_CG_OD2_diangle=geo.CA_CB_CG_OD2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    oxygen_d1= calculateCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle)
    OD1= Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    oxygen_d2= calculateCoordinates(CA, CB, CG, CG_OD2_length, CB_CG_OD2_angle, CA_CB_CG_OD2_diangle)
    OD2= Atom("OD2", oxygen_d2, 0.0, 1.0, " ", " OD2", 0, "O")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "ASP", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG) 
    res.add(OD1)
    res.add(OD2)
    return res

def makeAsn(segID,N, CA, C, O, geo):
    '''Creates an Asparagine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_OD1_length=geo.CG_OD1_length
    CB_CG_OD1_angle=geo.CB_CG_OD1_angle
    CA_CB_CG_OD1_diangle=geo.CA_CB_CG_OD1_diangle
    
    CG_ND2_length=geo.CG_ND2_length
    CB_CG_ND2_angle=geo.CB_CG_ND2_angle
    CA_CB_CG_ND2_diangle=geo.CA_CB_CG_ND2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    oxygen_d1= calculateCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle)
    OD1= Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    nitrogen_d2= calculateCoordinates(CA, CB, CG, CG_ND2_length, CB_CG_ND2_angle, CA_CB_CG_ND2_diangle)
    ND2= Atom("ND2", nitrogen_d2, 0.0, 1.0, " ", " ND2", 0, "N")
    res= Residue((' ', segID, ' '), "ASN", '    ')

    ##Create Residue Data Structure
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG) 
    res.add(OD1)
    res.add(ND2)
    return res

def makeGlu(segID, N, CA, C, O, geo):
    '''Creates a Glutamic Acid residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD_length=geo.CG_CD_length
    CB_CG_CD_angle=geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle=geo.CA_CB_CG_CD_diangle

    CD_OE1_length=geo.CD_OE1_length
    CG_CD_OE1_angle=geo.CG_CD_OE1_angle
    CB_CG_CD_OE1_diangle=geo.CB_CG_CD_OE1_diangle

    CD_OE2_length=geo.CD_OE2_length
    CG_CD_OE2_angle=geo.CG_CD_OE2_angle
    CB_CG_CD_OE2_diangle=geo.CB_CG_CD_OE2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle)
    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    oxygen_e1= calculateCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle)
    OE1= Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    oxygen_e2= calculateCoordinates(CB, CG, CD, CD_OE2_length, CG_CD_OE2_angle, CB_CG_CD_OE2_diangle)
    OE2= Atom("OE2", oxygen_e2, 0.0, 1.0, " ", " OE2", 0, "O")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "GLU", '    ')
    
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(OE1)
    res.add(OE2)
    return res

def makeGln(segID, N, CA, C, O, geo):
    '''Creates a Glutamine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD_length=geo.CG_CD_length
    CB_CG_CD_angle=geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle=geo.CA_CB_CG_CD_diangle
    
    CD_OE1_length=geo.CD_OE1_length
    CG_CD_OE1_angle=geo.CG_CD_OE1_angle
    CB_CG_CD_OE1_diangle=geo.CB_CG_CD_OE1_diangle
    
    CD_NE2_length=geo.CD_NE2_length
    CG_CD_NE2_angle=geo.CG_CD_NE2_angle
    CB_CG_CD_NE2_diangle=geo.CB_CG_CD_NE2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle)
    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    oxygen_e1= calculateCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle)
    OE1= Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    nitrogen_e2= calculateCoordinates(CB, CG, CD, CD_NE2_length, CG_CD_NE2_angle, CB_CG_CD_NE2_diangle)
    NE2= Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")


    ##Create Residue DS
    res= Residue((' ', segID, ' '), "GLN", '    ')
    
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(OE1)
    res.add(NE2)
    return res

def makeMet(segID, N, CA, C, O, geo):
    '''Creates a Methionine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_SD_length=geo.CG_SD_length
    CB_CG_SD_angle=geo.CB_CG_SD_angle
    CA_CB_CG_SD_diangle=geo.CA_CB_CG_SD_diangle
    
    SD_CE_length=geo.SD_CE_length
    CG_SD_CE_angle=geo.CG_SD_CE_angle
    CB_CG_SD_CE_diangle=geo.CB_CG_SD_CE_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    sulfur_d= calculateCoordinates(CA, CB, CG, CG_SD_length, CB_CG_SD_angle, CA_CB_CG_SD_diangle)
    SD= Atom("SD", sulfur_d, 0.0, 1.0, " ", " SD", 0, "S")
    carbon_e= calculateCoordinates(CB, CG, SD, SD_CE_length, CG_SD_CE_angle, CB_CG_SD_CE_diangle)
    CE= Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "MET", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(SD)
    res.add(CE)
    return res

def makeHis(segID, N, CA, C, O, geo):
    '''Creates a Histidine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_ND1_length=geo.CG_ND1_length
    CB_CG_ND1_angle=geo.CB_CG_ND1_angle
    CA_CB_CG_ND1_diangle=geo.CA_CB_CG_ND1_diangle
    
    CG_CD2_length=geo.CG_CD2_length
    CB_CG_CD2_angle=geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle=geo.CA_CB_CG_CD2_diangle
    
    ND1_CE1_length=geo.ND1_CE1_length
    CG_ND1_CE1_angle=geo.CG_ND1_CE1_angle
    CB_CG_ND1_CE1_diangle=geo.CB_CG_ND1_CE1_diangle
    
    CD2_NE2_length=geo.CD2_NE2_length
    CG_CD2_NE2_angle=geo.CG_CD2_NE2_angle
    CB_CG_CD2_NE2_diangle=geo.CB_CG_CD2_NE2_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    nitrogen_d1= calculateCoordinates(CA, CB, CG, CG_ND1_length, CB_CG_ND1_angle, CA_CB_CG_ND1_diangle)
    ND1= Atom("ND1", nitrogen_d1, 0.0, 1.0, " ", " ND1", 0, "N")
    carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle)
    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1= calculateCoordinates(CB, CG, ND1, ND1_CE1_length, CG_ND1_CE1_angle, CB_CG_ND1_CE1_diangle)
    CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    nitrogen_e2= calculateCoordinates(CB, CG, CD2, CD2_NE2_length, CG_CD2_NE2_angle, CB_CG_CD2_NE2_diangle)
    NE2= Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "HIS", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(ND1)
    res.add(CD2)
    res.add(CE1)
    res.add(NE2)
    return res

def makePro(segID, N, CA, C, O, geo):
    '''Creates a Proline residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle
    
    CG_CD_length=geo.CG_CD_length
    CB_CG_CD_angle=geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle=geo.CA_CB_CG_CD_diangle
    
    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d= calculateCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle)
    CD= Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")

    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), "PRO", '    ')
    
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)

    return res

def makePhe(segID, N, CA, C, O, geo):
    '''Creates a Phenylalanine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD1_length=geo.CG_CD1_length
    CB_CG_CD1_angle=geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle=geo.CA_CB_CG_CD1_diangle

    CG_CD2_length=geo.CG_CD2_length
    CB_CG_CD2_angle=geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle= geo.CA_CB_CG_CD2_diangle
    
    CD1_CE1_length=geo.CD1_CE1_length
    CG_CD1_CE1_angle=geo.CG_CD1_CE1_angle
    CB_CG_CD1_CE1_diangle=geo.CB_CG_CD1_CE1_diangle

    CD2_CE2_length=geo.CD2_CE2_length
    CG_CD2_CE2_angle=geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle=geo.CB_CG_CD2_CE2_diangle

    CE1_CZ_length=geo.CE1_CZ_length
    CD1_CE1_CZ_angle=geo.CD1_CE1_CZ_angle
    CG_CD1_CE1_CZ_diangle=geo.CG_CD1_CE1_CZ_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle)
    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle)
    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1= calculateCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle)
    CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle)
    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_z= calculateCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle)
    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")

    ##Create Residue Data Structures
    res= Residue((' ', segID, ' '), "PHE", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CE1)
    res.add(CD2)
    res.add(CE2)
    res.add(CZ)
    return res

def makeTyr(segID, N, CA, C, O, geo):
    '''Creates a Tyrosine residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle
    
    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD1_length=geo.CG_CD1_length
    CB_CG_CD1_angle=geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle=geo.CA_CB_CG_CD1_diangle
    
    CG_CD2_length=geo.CG_CD2_length
    CB_CG_CD2_angle=geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle=geo.CA_CB_CG_CD2_diangle
    
    CD1_CE1_length=geo.CD1_CE1_length
    CG_CD1_CE1_angle=geo.CG_CD1_CE1_angle
    CB_CG_CD1_CE1_diangle=geo.CB_CG_CD1_CE1_diangle

    CD2_CE2_length=geo.CD2_CE2_length
    CG_CD2_CE2_angle=geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle=geo.CB_CG_CD2_CE2_diangle

    CE1_CZ_length=geo.CE1_CZ_length
    CD1_CE1_CZ_angle=geo.CD1_CE1_CZ_angle
    CG_CD1_CE1_CZ_diangle=geo.CG_CD1_CE1_CZ_diangle

    CZ_OH_length=geo.CZ_OH_length
    CE1_CZ_OH_angle=geo.CE1_CZ_OH_angle
    CD1_CE1_CZ_OH_diangle=geo.CD1_CE1_CZ_OH_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle)
    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle)
    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1= calculateCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle)
    CE1= Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle)
    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_z= calculateCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle)
    CZ= Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    oxygen_h= calculateCoordinates(CD1, CE1, CZ, CZ_OH_length, CE1_CZ_OH_angle, CD1_CE1_CZ_OH_diangle)
    OH= Atom("OH", oxygen_h, 0.0, 1.0, " ", " OH", 0, "O")

    ##Create Residue Data S
    res= Residue((' ', segID, ' '), "TYR", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CE1)
    res.add(CD2)
    res.add(CE2)
    res.add(CZ)
    res.add(OH)
    return res

def makeTrp(segID, N, CA, C, O, geo):
    '''Creates a Tryptophan residue'''
    ##R-Group
    CA_CB_length=geo.CA_CB_length
    C_CA_CB_angle=geo.C_CA_CB_angle
    N_C_CA_CB_diangle=geo.N_C_CA_CB_diangle

    CB_CG_length=geo.CB_CG_length
    CA_CB_CG_angle=geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle=geo.N_CA_CB_CG_diangle

    CG_CD1_length=geo.CG_CD1_length
    CB_CG_CD1_angle=geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle=geo.CA_CB_CG_CD1_diangle

    CG_CD2_length=geo.CG_CD2_length
    CB_CG_CD2_angle=geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle=geo.CA_CB_CG_CD2_diangle
    
    CD1_NE1_length=geo.CD1_NE1_length
    CG_CD1_NE1_angle=geo.CG_CD1_NE1_angle
    CB_CG_CD1_NE1_diangle=geo.CB_CG_CD1_NE1_diangle

    CD2_CE2_length=geo.CD2_CE2_length
    CG_CD2_CE2_angle=geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle=geo.CB_CG_CD2_CE2_diangle

    CD2_CE3_length=geo.CD2_CE3_length
    CG_CD2_CE3_angle=geo.CG_CD2_CE3_angle
    CB_CG_CD2_CE3_diangle=geo.CB_CG_CD2_CE3_diangle

    CE2_CZ2_length=geo.CE2_CZ2_length
    CD2_CE2_CZ2_angle=geo.CD2_CE2_CZ2_angle
    CG_CD2_CE2_CZ2_diangle=geo.CG_CD2_CE2_CZ2_diangle

    CE3_CZ3_length=geo.CE3_CZ3_length
    CD2_CE3_CZ3_angle=geo.CD2_CE3_CZ3_angle
    CG_CD2_CE3_CZ3_diangle=geo.CG_CD2_CE3_CZ3_diangle

    CZ2_CH2_length=geo.CZ2_CH2_length
    CE2_CZ2_CH2_angle=geo.CE2_CZ2_CH2_angle
    CD2_CE2_CZ2_CH2_diangle=geo.CD2_CE2_CZ2_CH2_diangle

    carbon_b= calculateCoordinates(N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle)
    CB= Atom("CB", carbon_b, 0.0 , 1.0, " "," CB", 0,"C")
    carbon_g= calculateCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle)
    CG= Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1= calculateCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle)
    CD1= Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2= calculateCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle)
    CD2= Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    nitrogen_e1= calculateCoordinates(CB, CG, CD1, CD1_NE1_length, CG_CD1_NE1_angle, CB_CG_CD1_NE1_diangle)
    NE1= Atom("NE1", nitrogen_e1, 0.0, 1.0, " ", " NE1", 0, "N")
    carbon_e2= calculateCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle)
    CE2= Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_e3= calculateCoordinates(CB, CG, CD2, CD2_CE3_length, CG_CD2_CE3_angle, CB_CG_CD2_CE3_diangle)
    CE3= Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")

    carbon_z2= calculateCoordinates(CG, CD2, CE2, CE2_CZ2_length, CD2_CE2_CZ2_angle, CG_CD2_CE2_CZ2_diangle)
    CZ2= Atom("CZ2", carbon_z2, 0.0, 1.0, " ", " CZ2", 0, "C")

    carbon_z3= calculateCoordinates(CG, CD2, CE3, CE3_CZ3_length, CD2_CE3_CZ3_angle, CG_CD2_CE3_CZ3_diangle)
    CZ3= Atom("CZ3", carbon_z3, 0.0, 1.0, " ", " CZ3", 0, "C")

    carbon_h2= calculateCoordinates(CD2, CE2, CZ2, CZ2_CH2_length, CE2_CZ2_CH2_angle, CD2_CE2_CZ2_CH2_diangle)
    CH2= Atom("CH2", carbon_h2, 0.0, 1.0, " ", " CH2", 0, "C")
    
    ##Create Residue DS
    res= Residue((' ', segID, ' '), "TRP", '    ')
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CD2)

    res.add(NE1)
    res.add(CE2)
    res.add(CE3)

    res.add(CZ2)
    res.add(CZ3)

    res.add(CH2)
    return res

def initialize_res(residue):
    '''Creates a new structure containing a single amino acid. The type and
    geometry of the amino acid are determined by the argument, which has to be
    either a geometry object or a single-letter amino acid code.
    The amino acid will be placed into chain A of model 0.'''
    
    if isinstance( residue, Geo ):
        geo = residue
    else:
        geo=geometry(residue) 
    
    segID=1
    AA= geo.residue_name
    CA_N_length=geo.CA_N_length
    CA_C_length=geo.CA_C_length
    N_CA_C_angle=geo.N_CA_C_angle
    
    CA_coord= numpy.array([0.,0.,0.])
    C_coord= numpy.array([CA_C_length,0,0])
    N_coord = numpy.array([CA_N_length*math.cos(N_CA_C_angle*(math.pi/180.0)),CA_N_length*math.sin(N_CA_C_angle*(math.pi/180.0)),0])

    N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")
    CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")
    C= Atom("C", C_coord, 0.0, 1.0, " ", " C",0,"C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length=geo.C_O_length
    CA_C_O_angle=geo.CA_C_O_angle
    N_CA_C_O_diangle=geo.N_CA_C_O_diangle
    
    carbonyl=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle)
    O= Atom("O",carbonyl , 0.0 , 1.0, " "," O", 0, "O")

    if(AA=='G'):
        res=makeGly(segID, N, CA, C, O, geo)
    elif(AA=='A'):
        res=makeAla(segID, N, CA, C, O, geo)
    elif(AA=='S'):
        res=makeSer(segID, N, CA, C, O, geo)
    elif(AA=='C'):
        res=makeCys(segID, N, CA, C, O, geo)
    elif(AA=='V'):
        res=makeVal(segID, N, CA, C, O, geo)
    elif(AA=='I'):
        res=makeIle(segID, N, CA, C, O, geo)
    elif(AA=='L'):
        res=makeLeu(segID, N, CA, C, O, geo)
    elif(AA=='T'):
        res=makeThr(segID, N, CA, C, O, geo)
    elif(AA=='R'):
        res=makeArg(segID, N, CA, C, O, geo)
    elif(AA=='K'):
        res=makeLys(segID, N, CA, C, O, geo)
    elif(AA=='D'):
        res=makeAsp(segID, N, CA, C, O, geo)
    elif(AA=='E'):
        res=makeGlu(segID, N, CA, C, O, geo)
    elif(AA=='N'):
        res=makeAsn(segID, N, CA, C, O, geo)
    elif(AA=='Q'):
        res=makeGln(segID, N, CA, C, O, geo)
    elif(AA=='M'):
        res=makeMet(segID, N, CA, C, O, geo)
    elif(AA=='H'):
        res=makeHis(segID, N, CA, C, O, geo)
    elif(AA=='P'):
        res=makePro(segID, N, CA, C, O, geo)
    elif(AA=='F'):
        res=makePhe(segID, N, CA, C, O, geo)
    elif(AA=='Y'):
        res=makeTyr(segID, N, CA, C, O, geo)
    elif(AA=='W'):
        res=makeTrp(segID, N, CA, C, O, geo)
    else:
        res=makeGly(segID, N, CA, C, O, geo)

    cha= Chain('A')
    cha.add(res)
    
    mod= Model(0)
    mod.add(cha)

    struc= Structure('X')
    struc.add(mod)
    return struc


def getReferenceResidue(structure):
    '''Returns the last residue of chain A model 0 of the given structure.
    
    This function is a helper function that should not normally be called
    directly.'''

    # If the following line doesn't work we're in trouble.
    # Likely initialize_res() wasn't called.
    resRef = structure[0]['A'].child_list[-1]
    
    # If the residue is not an amino acid we're in trouble.
    # Likely somebody is trying to append residues to an existing
    # structure that has non-amino-acid molecules in the chain.
    assert is_aa(resRef)
        
    return resRef

def add_residue_from_geo(structure, geo):
    '''Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added is determined by
    the geometry object given as second argument.
    
    This function is a helper function and should not normally be called
    directly. Call add_residue() instead.'''
    resRef= getReferenceResidue(structure)
    AA=geo.residue_name
    segID= resRef.get_id()[1]
    segID+=1

    ##geometry to bring together residue
    peptide_bond=geo.peptide_bond
    CA_C_N_angle=geo.CA_C_N_angle
    C_N_CA_angle=geo.C_N_CA_angle

    ##Backbone Coordinages
    N_CA_C_angle=geo.N_CA_C_angle
    CA_N_length=geo.CA_N_length
    CA_C_length=geo.CA_C_length
    phi= geo.phi
    psi_im1=geo.psi_im1
    omega=geo.omega

    N_coord=calculateCoordinates(resRef['N'], resRef['CA'], resRef['C'], peptide_bond, CA_C_N_angle, psi_im1)
    N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")

    CA_coord=calculateCoordinates(resRef['CA'], resRef['C'], N, CA_N_length, C_N_CA_angle, omega)
    CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")

    C_coord=calculateCoordinates(resRef['C'], N, CA, CA_C_length, N_CA_C_angle, phi)
    C= Atom("C", C_coord, 0.0, 1.0, " ", " C",0,"C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length=geo.C_O_length
    CA_C_O_angle=geo.CA_C_O_angle
    N_CA_C_O_diangle=geo.N_CA_C_O_diangle

    carbonyl=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle)
    O= Atom("O",carbonyl , 0.0 , 1.0, " "," O", 0, "O")
    
    if(AA=='G'):
        res=makeGly(segID, N, CA, C, O, geo)
    elif(AA=='A'):
        res=makeAla(segID, N, CA, C, O, geo)
    elif(AA=='S'):
        res=makeSer(segID, N, CA, C, O, geo)
    elif(AA=='C'):
        res=makeCys(segID, N, CA, C, O, geo)
    elif(AA=='V'):
        res=makeVal(segID, N, CA, C, O, geo)
    elif(AA=='I'):
        res=makeIle(segID, N, CA, C, O, geo)
    elif(AA=='L'):
        res=makeLeu(segID, N, CA, C, O, geo)
    elif(AA=='T'):
        res=makeThr(segID, N, CA, C, O, geo)
    elif(AA=='R'):
        res=makeArg(segID, N, CA, C, O, geo)
    elif(AA=='K'):
        res=makeLys(segID, N, CA, C, O, geo)
    elif(AA=='D'):
        res=makeAsp(segID, N, CA, C, O, geo)
    elif(AA=='E'):
        res=makeGlu(segID, N, CA, C, O, geo)
    elif(AA=='N'):
        res=makeAsn(segID, N, CA, C, O, geo)
    elif(AA=='Q'):
        res=makeGln(segID, N, CA, C, O, geo)
    elif(AA=='M'):
        res=makeMet(segID, N, CA, C, O, geo)
    elif(AA=='H'):
        res=makeHis(segID, N, CA, C, O, geo)
    elif(AA=='P'):
        res=makePro(segID, N, CA, C, O, geo)
    elif(AA=='F'):
        res=makePhe(segID, N, CA, C, O, geo)
    elif(AA=='Y'):
        res=makeTyr(segID, N, CA, C, O, geo)
    elif(AA=='W'):
        res=makeTrp(segID, N, CA, C, O, geo)
    else:
        res=makeGly(segID, N, CA, C, O, geo)
        
    resRef['O'].set_coord(calculateCoordinates(res['N'], resRef['CA'], resRef['C'], C_O_length, CA_C_O_angle, 180.0))

    ghost= Atom("N", calculateCoordinates(res['N'], res['CA'], res['C'], peptide_bond, CA_C_N_angle, psi_im1), 0.0 , 0.0, " ","N", 0, "N")
    res['O'].set_coord(calculateCoordinates( res['N'], res['CA'], res['C'], C_O_length, CA_C_O_angle, 180.0))

    structure[0]['A'].add(res)
    return structure


def make_extended_structure(AA_chain):
    '''Place a sequence of amino acids into a peptide in the extended
    conformation. The argument AA_chain holds the sequence of amino
    acids to be used.'''
    geo = geometry(AA_chain[0])
    struc=initialize_res(geo)
    
    for i in range(1,len(AA_chain)): 
        AA = AA_chain[i]
        geo = geometry(AA)
        add_residue(struc, geo)

    return struc

    
def add_residue(structure, residue, phi=-120, psi_im1=140, omega=-370):
    '''Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added can be specified
    in two ways: either as a geometry object (in which case
    the remaining arguments phi, psi_im1, and omega are ignored) or as a
    single-letter amino-acid code. In the latter case, the optional
    arguments phi, psi_im1, and omega specify the corresponding backbone
    angles.
    
    When omega is specified, it needs to be a value greater than or equal
    to -360. Values below -360 are ignored.''' 
    
    if isinstance( residue, Geo ):
        geo = residue
    else:
        geo=geometry(residue) 
        geo.phi=phi
        geo.psi_im1=psi_im1
        if omega>-361:
            geo.omega=omega
    
    add_residue_from_geo(structure, geo)
    return structure
    

    
def make_structure(AA_chain,phi,psi_im1,omega=[]):
    '''Place a sequence of amino acids into a peptide with specified
    backbone dihedral angles. The argument AA_chain holds the
    sequence of amino acids to be used. The arguments phi and psi_im1 hold
    lists of backbone angles, one for each amino acid, *starting from
    the second amino acid in the chain*. The argument 
    omega (optional) holds a list of omega angles, also starting from
    the second amino acid in the chain.'''
    geo = geometry(AA_chain[0])
    struc=initialize_res(geo)

    if len(omega)==0:
        for i in range(1,len(AA_chain)): 
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i-1], psi_im1[i-1])
    else:
        for i in range(1,len(AA_chain)): 
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i-1], psi_im1[i-1], omega[i-1])

    return struc
    
    
def make_structure_from_geos(geos):
    '''Creates a structure out of a list of geometry objects.'''
    model_structure=initialize_res(geos[0])
    for i in range(1,len(geos)):
        model_structure=add_residue(model_structure, geos[i])

    return model_structure

