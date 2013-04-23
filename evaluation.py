#!/usr/bin/python

from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Structure import *
from Bio.PDB.Vector import *
from Bio.PDB.Entity import*
import math
import Geometry
import PeptideBuilder
import numpy

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
	    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
	    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
	    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

def build_linear_model(pdb_filename):
    parser=PDBParser()
    structure=parser.get_structure( 'sample', pdb_filename )
    model=structure[0]
    chain=model['A']
    model_structure_geo=[]
    for res in chain:
        if(res.get_resname() in resdict.keys()):
            tempgeo=Geometry.geometry(resdict[res.get_resname()])
            model_structure_geo.append(tempgeo)
    model_structure=PeptideBuilder.initialize_res(model_structure_geo[0])
    for i in range(1,len(model_structure_geo)):
        model_structure=PeptideBuilder.add_residue(model_structure, model_structure_geo[i])

    return model_structure                

def make_pdb_file(struct, file_nom):
    outfile = PDBIO()
    outfile.set_structure(struct)
    outfile.save( file_nom)
    return file_nom        

def build_backbone_model(pdb_filename):
    parser=PDBParser()
    structure=parser.get_structure( 'sample', pdb_filename )
    model=structure[0]
    chain=model['A']
    model_structure_geo=[]
    prev="0"
    N_prev="0"
    CA_prev="0"
    CO_prev="0"
    ##O_prev="0"
    prev_res=""
    rad=180.0/math.pi
    for res in chain:
        if(res.get_resname() in resdict.keys()):
            geo=Geometry.geometry(resdict[res.get_resname()])
            if(prev=="0"):
                 N_prev=res['N']
                 CA_prev=res['CA']
                 C_prev=res['C']
                 ##O_prev=res['O']
                 prev="1"
            else:
                 n1=N_prev.get_vector()
                 ca1=CA_prev.get_vector()
                 c1=C_prev.get_vector()
                 ##o1=O_prev.get_vector()
                            
                 ##O_curr=res['O']
                 C_curr=res['C']
                 N_curr=res['N']
                 CA_curr=res['CA']
                                             
                 ##o=O_curr.get_vector()
                 c=C_curr.get_vector()
                 n=N_curr.get_vector()
                 ca=CA_curr.get_vector()

                 geo.CA_C_N_angle=calc_angle(ca1, c1, n)*rad
                 geo.C_N_CA_angle=calc_angle(c1, n, ca)*rad
                 geo.CA_N_length= CA_curr-N_curr
                 geo.CA_C_length= CA_curr-C_curr
                 geo.peptide_bond= N_curr-C_prev

                 psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                 omega= calc_dihedral(ca1, c1, n, ca) ##goes to current res
                 phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                 geo.psi=psi*rad
                 geo.omega=omega*rad
                 geo.phi=phi*rad

                 geo.CA_N_length= CA_curr - N_curr
                 geo.CA_C_length= CA_curr - C_curr 
                 ##geo.C_O_length= C_curr - O_curr

                 geo.N_CA_C_angle= calc_angle(n, ca, c)*rad
                 ##geo.CA_C_O_angle= calc_angle(ca, c, o)*rad

                 ##geo.N_CA_C_O= calc_dihedral(n, ca, c, o)*rad

                 N_prev=res['N']
                 CA_prev=res['CA']
                 C_prev=res['C']
                 ##O_prev=res['O']
                               
                        
            model_structure_geo.append(geo)
##    model_structure=PeptideBuilder.initialize_res(model_structure_geo[0])
##    for i in range(1,len(model_structure_geo)):
##        model_structure=PeptideBuilder.add_residue(model_structure, model_structure_geo[i])
##
##    return model_structure
    return model_structure_geo

def build_limited_backbone_model(pdb_filename):
    parser=PDBParser()
    structure=parser.get_structure( 'sample', pdb_filename )
    model=structure[0]
    chain=model['A']
    model_structure_geo=[]
    prev="0"
    N_prev="0"
    CA_prev="0"
    CO_prev="0"
    prev_res=""
    rad=180.0/math.pi
    for res in chain:
        if(res.get_resname() in resdict.keys()):
            geo=Geometry.geometry(resdict[res.get_resname()])
            if(prev=="0"):
                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']
                prev="1"
            else:
                n1=N_prev.get_vector()
                ca1=CA_prev.get_vector()
                c1=C_prev.get_vector()
                                
                C_curr=res['C']
                N_curr=res['N']
                CA_curr=res['CA']
                                                
                c=C_curr.get_vector()
                n=N_curr.get_vector()
                ca=CA_curr.get_vector()

                geo.CA_C_N_angle=calc_angle(ca1, c1, n)*rad
                geo.C_N_CA_angle=calc_angle(c1, n, ca)*rad

                psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                omega= calc_dihedral(ca1, c1, n, ca) ##goes to current res
                phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                geo.psi=psi*rad
                geo.omega=omega*rad
                geo.phi=phi*rad

                geo.N_CA_C_angle= calc_angle(n, ca, c)*rad
                ##geo.CA_C_O_angle= calc_angle(ca, c, o)*rad

                ##geo.N_CA_C_O= calc_dihedral(n, ca, c, o)*rad

                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']
                ##O_prev=res['O']
                                
                        
            model_structure_geo.append(geo)
##        model_structure=PeptideBuilder.initialize_res(model_structure_geo[0])
##        for i in range(1,len(model_structure_geo)):
##                model_structure=PeptideBuilder.add_residue(model_structure, model_structure_geo[i])
##
##        return model_structure
    return model_structure_geo

        
def build_rough_model(pdb_filename):
    parser=PDBParser()
    structure=parser.get_structure( 'sample', pdb_filename )
    model=structure[0]
    chain=model['A']
    seq=""
    phi_diangle=[]
    psi_diangle=[]
    omega_diangle=[]
    for res in chain:
        if(res.get_resname() in resdict.keys()):
                        
            seq+=resdict[res.get_resname()]
            if(len(seq)==1):
                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']
            else:   
                n1=N_prev.get_vector()
                ca1=CA_prev.get_vector()
                c1=C_prev.get_vector()
                               
                C_curr=res['C']
                N_curr=res['N']
                CA_curr=res['CA']
                                                
                c=C_curr.get_vector()
                n=N_curr.get_vector()
                ca=CA_curr.get_vector()

                psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                omega= calc_dihedral(ca1, c1, n, ca)
                phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                phi_diangle.append(phi*180.0/math.pi)
                psi_diangle.append(psi*180.0/math.pi)
                omega_diangle.append(omega*180.0/math.pi)

                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']
                                
    model_structure_omega= PeptideBuilder.make_structure( seq, phi_diangle, psi_diangle, omega_diangle )
    model_structure_rough= PeptideBuilder.make_structure( seq, phi_diangle, psi_diangle)
    return model_structure_omega, model_structure_rough   
##geo = Geometry.geometry("G")
##structure = PeptideBuilder.initialize_res(geo)
##for i in range(4):
##	structure = PeptideBuilder.add_residue(structure, geo)

def compare_structure( reference, alternate):
    parser=PDBParser()

    ref_struct=parser.get_structure('Reference', reference)
    alt_struct= parser.get_structure("Alternate", alternate)


    ref_model=ref_struct[0]
    ref_chain=ref_model['A']

    alt_model=alt_struct[0]
    alt_chain=alt_model['A']

    ref_atoms=[]
    alt_atoms=[]

    for ref_res in ref_chain:
        if(ref_res.get_resname() in resdict.keys()):
            ref_atoms.append(ref_res['CA'])

    for alt_res in alt_chain:
        if(alt_res.get_resname() in resdict.keys()):
             alt_atoms.append(alt_res['CA'])

    super_imposer= Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(alt_model.get_atoms())

    make_pdb_file(alt_struct, "Aligned_" + alternate)

    full= super_imposer.rms

    super_imposer_50= Superimposer()
    super_imposer_50.set_atoms(ref_atoms[:50], alt_atoms[:50])
    super_imposer_50.apply(alt_model.get_atoms())

    make_pdb_file(alt_struct, "Aligned_50_" + alternate)

    f_50= super_imposer_50.rms

    super_imposer_150= Superimposer()
    super_imposer_150.set_atoms(ref_atoms[:150], alt_atoms[:150])
    super_imposer_150.apply(alt_model.get_atoms())

    make_pdb_file(alt_struct, "Aligned_150_" + alternate)

    f_150= super_imposer_150.rms

    return f_50, f_150, full, len(ref_atoms)

        
        
def test_PeptideBuilder(pdb_File):
    ##structure_backbone= PeptideBuilder.make_geo_structure(build_backbone_model(pdb_File))
    ##structure_omega, structure_rough=build_rough_model(pdb_File)
    ##structure_limited= PeptideBuilder.make_geo_structure(build_limited_backbone_model(pdb_File))

    ##RMS_backbone_50,RMS_backbone_150, RMS_backbone, size= compare_structure(pdb_File, make_pdb_file(structure_backbone, "Backbone_" + pdb_File))
    ##RMS_rough_50, RMS_rough_150, RMS_rough, size= compare_structure(pdb_File, make_pdb_file(structure_rough, "Rough_" + pdb_File))
    ##RMS_omega_50, RMS_omega_150, RMS_omega, size= compare_structure(pdb_File, make_pdb_file(structure_omega, "Omega_" + pdb_File))
    RMS_limited_50, RMS_limited_150, RMS_limited, size= compare_structure(pdb_File, make_pdb_file(structure_limited, "Limited_" + pdb_File))
    ##output_line= "%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%i\n" % (pdb_File[:-4],RMS_rough_50, RMS_omega_50, RMS_backbone_50,RMS_rough_150, RMS_omega_150, RMS_backbone_150, RMS_rough, RMS_omega, RMS_backbone, size)
    output_line= "%s\t%0.2f\t%0.2f\t%0.2f\t%i\n" % (pdb_File[:-4],RMS_limited_50, RMS_limited_150, RMS_limited, size)
    return output_line
##"1LMB.pdb"
##"2jfg.pdb"
##"3cmo.pdb"
test_files=["4eoi.pdb","3cuq.pdb","1lnm.pdb","2ocf.pdb","2o6r.pdb","1nbw.pdb","7tim.pdb", "1vca.pdb", "1aq7.pdb", "1gfl.pdb"]
##<100 "2b88.pdb","1ayi.pdb","2faz.pdb","3rcp.pdb","1mi0.pdb", "1fme.pdb", "2jof.pdb",
##<150 "2jab.pdb","1ves.pdb","3q0h.pdb","2etd.pdb","2lyz.pdb","3sdh.pdb", "3ecf.pdb",
f_out=open("Evaluation_PeptideBuilder_limited.txt","w")
##f_out.write("PDB\tPhi-Psi-50\tPhi-Psi-Omega-50\tAll_backbone-50\tPhi-Psi-150\tPhi-Psi-Omega-150\tAll_backbone-150\tPhi-Psi\tPhi-Psi-Omega\tAll_backbone\tlength\n")
f_out.write("PDB\tlimited50\tlimited150\tlimited_full\tlength\n")
for i in test_files:
    print i
    f_out.write(test_PeptideBuilder(i))

f_out.close()
