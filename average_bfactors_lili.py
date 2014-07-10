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

for residue in model["C"] :
    average_bfactors[ residue.get_id()[1] ] = 0.0
    
for chain in model.get_list():
    for residue in chain.get_list():
      if residue.has_id("CA"):
        ca=residue["CA"]
        average_bfactors[ residue.get_id()[1] ] += float( ca.get_bfactor() ) / float( len( model.get_list() ) )

for chain in model.get_list():
    for residue in chain.get_list():
        for atom in residue.get_list():
            atom.set_bfactor( average_bfactors[ residue.get_id()[1] ] )

w = PDBIO()
w.set_structure(structure)
w.save('sym_'+options.pdbfile)

