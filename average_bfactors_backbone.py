#!/usr/bin/python

# Xavier Ambroggio, 2012
# ambroggiox@niaid.nih.gov

from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-f", "--pdbfile", default = None, type = "string", help="pdb structure file for additional 3-coord cartesian per residue")
(options,args)=parser.parse_args()

parser=PDBParser()

structure=parser.get_structure("mystruct", options.pdbfile)
model=structure[0]

average_bfactors = {}

for residue in model["X"] :
    average_bfactors[ residue.get_id()[1] ] = 0.0
    
for chain in model.get_list():
    for residue in chain.get_list():
      sum_over_natoms = 0  
      if residue.has_id("N"):
        ca=residue["N"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("H"):
        ca=residue["H"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("H1"):
        ca=residue["H1"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("CA"):
        ca=residue["CA"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("C"):
        ca=residue["C"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("O"):
        ca=residue["O"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("OXT"):
        ca=residue["OXT"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) 
        sum_over_natoms += 1
      if residue.has_id("CA"):
          average_bfactors[  residue.get_id()[1] ] /= sum_over_natoms

for chain in model.get_list():
    for residue in chain.get_list():
        for atom in residue.get_list():
            atom.set_bfactor( average_bfactors[ residue.get_id()[1] ] )

w = PDBIO()
w.set_structure(structure)
w.save('avg_b_backbone_'+options.pdbfile)

