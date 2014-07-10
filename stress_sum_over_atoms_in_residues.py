#!/bin/python
from numpy import *
#from pylab import *
import sys
import re
from optparse import OptionParser
import os
import gzip
from scipy.stats import gaussian_kde
#import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser



if __name__ == "__main__":
   parser=OptionParser()
   parser.add_option("-s", "--structure", default=None, type="string", help="pdb file")
   parser.add_option("-f", "--file", default=None, type="string", help="stress .dat file to process")
   parser.add_option("-S", "--smooth", default=1, type="int", help="number of consecutive points to use for smoothing")
   parser.add_option("-n", "--num_blocks", default=1, type="int",help="number of subsets to split simulation into")
   (options,args)=parser.parse_args()

   #use biopython to iterate over residues, grab atom numbers in residue, sum up stresses, output to xvg file

   reslistfile = open(options.structure + ".reslist","w")
   stress_datafile = open(options.file, "r")
   print "loading data"
   numlines = 0
   
   p = PDBParser()
   structure = p.get_structure('X', options.structure)
   num_atoms = 0
   print "loaded PDB file"
   for model in structure:
       for chain in model:
           for residue in chain:
              for atom in residue:
                num_atoms += 1

   print "number of atoms: "+str(num_atoms)
           
   for line in stress_datafile:
      if (numlines == 0):
         numfields=len(line.split())
      numlines += 1
      #print "numlines :"+str(numlines)
      if (numlines % 10000 == 0):
         print "line of stress .dat file: "+str(numlines)
   print "line of stress .dat file: "+str(numlines)
   
   stress_data = zeros( (numlines / options.smooth, numfields), float32)
   stress_datafile.seek(0) #rewind
   
   i = 0
   j = 1
   for line in stress_datafile:
      if(i > 0 and i % j == 0): 
         stress_data[(i - 1) / options.smooth, :] /= j #average
         j = 1 #reset counter
      thisline = zeros((numfields),float64)
      thisline[:] = line.split()
      stress_data[i / options.smooth,:] += thisline 
      i += 1
      j += 1
      if (i % 10000 == 0):
         print "line of stress .dat file: "+str(i)
   print "line of stress .dat file: "+str(i)
   
   print "finished loading data"
   stress_xvg_dir = str(options.file)+"_blocks_"+str(options.num_blocks)+"_smooth_"+str(options.smooth)+"_stress_xvg"
                              
   number_of_stress_lines = stress_data.shape[0]
   number_of_stress_lines_per_block = int(number_of_stress_lines / options.num_blocks)
   stress_file_sub_blocks_right_breaks = zeros((options.num_blocks), int64)
   stress_file_sub_blocks_left_breaks = zeros((options.num_blocks), int64)
   
   right_break = 0
   os.system("mkdir "+stress_xvg_dir)
   for i in range(options.num_blocks):
       os.system("mkdir "+stress_xvg_dir+"/run"+str(i+1))
       stress_file_sub_blocks_left_breaks[i] = right_break
       right_break +=  number_of_stress_lines_per_block
       stress_file_sub_blocks_right_breaks[i] = right_break

   print "left breaks:"
   print  stress_file_sub_blocks_left_breaks
   print "right breaks:"
   print  stress_file_sub_blocks_right_breaks

   
   
   print "loaded PDB file"
   for model in structure:
       ires = 0
       for chain in model:
           for residue in chain:
               ires += 1
               stress_this_residue = zeros((options.num_blocks, number_of_stress_lines_per_block), float64)
               resname = residue.get_resname()
               resnum = int((residue.get_full_id()[-1])[1])
               chain = residue.get_full_id()[-2]
               reslistfile.write(str(resnum)+" "+str(resname)+" "+str(ires)+str(chain)+"\n")
               xvgfiles = []
               print "residue: "+str(resnum)+" "+str(resname)
               for block in range(options.num_blocks):
                   xvgfiles.append( open(str(stress_xvg_dir)+"/run"+str(block+1)+"/chi1"+str(resname)+str(ires)+".xvg","w"))
               for atom in residue:
                   atom_num = atom.get_serial_number()
                   for block in range(options.num_blocks):                       
                       stress_this_residue[block,:] += stress_data[stress_file_sub_blocks_left_breaks[block]:stress_file_sub_blocks_right_breaks[block], atom_num - 1]  #  performs sum over atoms in this residue as iterations progress, atom_num - 1 to avoid off-by-one error
               for block in range(options.num_blocks):
                   for i in range( number_of_stress_lines_per_block ):
                       xvgfiles[block].write(str(i)+" "+str(stress_this_residue[block,i])+"\n")
                   xvgfiles[block].close()




                       



