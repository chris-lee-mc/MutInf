#!/usr/bin/python

from dihedral_mutent import *
from PDBlite import *

def compute_CA_dist_matrix(reslist,pdb_obj):  #pdb_obj is from PDBlite.py                                                                                     
    #print "reslist:"
    #print reslist
    CA_dist_matrix = zeros((len(reslist),len(reslist)),float32)
    #print "CA dist matrix:"
    #print CA_dist_matrix
    for res_ind1, myres1 in zip(range(len(reslist)), reslist):
        for res_ind2, myres2 in zip(range(res_ind1, len(reslist)), reslist[res_ind1:]):
            #print "\n#### Working on residues %s and %s (%s and %s):" % (myres1.num, myres2.num, myres1.name, myres2.name)
            BBi = pdb_obj.get(myres1.chain,myres1.num)
            BBj = pdb_obj.get(myres2.chain,myres2.num)
            CAi = BBi.get_atom("CA") #CA alpha carbon from res i                                                                                           
            CAj = BBj.get_atom("CA") #CA alpha carbon from res j                                                                                           
            CA_dist_matrix[res_ind1,res_ind2] = sqrt(CAi.calc_dist2(CAj))
            CA_dist_matrix[res_ind2,res_ind1] = CA_dist_matrix[res_ind1,res_ind2] # for symmetry                                                           

    return CA_dist_matrix


#main

try:
   import run_profile
   run_profile.fix_args()
except: pass

usage="%prog [-s mypdb.pdb] [--matrix mystuff_mutinf_res_sum_0diag.txt] resfile  # where resfile is in the format <1-based-index> <aa type> <res num>"
parser=OptionParser(usage)
parser.add_option("-s", "--structure", default=None, type="string", help="filenames to load PDB from")
parser.add_option("-m", "--matrix",    default=None, type="string", help="res x res mutinf matrix to read from")
parser.add_option("-d", "--distance",  default=5.0, type="float", help="cutoff distance for filtering MutInf matrix")
(options,args)=parser.parse_args()

if len(filter(lambda x: x==None, (options.structure))) == 1:
    parser.error("ERROR exactly one of --structure must be specified")

resfile_fn = args[0]

print options
print "resfile fn: "+resfile_fn

print "structure:"+str(options.structure)

run_params = RunParameters(num_structs=1,resfile_fn=resfile_fn)
pdb_lines = commands.getoutput("egrep '^ATOM ' %s" % options.structure).split("\n")
#print pdb_lines
pdb_obj = PDBlite.PDB(pdb_lines,fn=options.structure)
print pdb_obj
reslist = load_resfile(run_params,load_angles=False)
#print "reslist:"
#print reslist
CA_dist_matrix = compute_CA_dist_matrix(reslist,pdb_obj)
print CA_dist_matrix
print CA_dist_matrix.shape
name_num_list=[]
for res in reslist: name_num_list.append(res.name + str(res.num))
output_matrix(resfile_fn+"_CA_dist_matrix.txt",CA_dist_matrix,name_num_list,name_num_list)

if(options.matrix != None):
    mutinf_res_matrix, rownames, colnames = read_res_matrix(options.matrix)
    # now, filter res matrix at various fractions of kT
    # start with (S/k * kT) >= 1/2 kT
    print mutinf_res_matrix.shape
    half_kT_mutinfs = copy(mutinf_res_matrix)
    half_kT_mutinfs[half_kT_mutinfs < 0.5 * (8.314*300/(4.184*1000)) ] = 0
    #next, take (S/k * kT) >= 1/4 kT
    quarter_kT_mutinfs = array(mutinf_res_matrix)
    quarter_kT_mutinfs[quarter_kT_mutinfs < 0.25 * (8.314*300/(4.184*1000)) ] = 0
    #next, take (S/k * kT) >= 0.1 kT
    tenth_kT_mutinfs = array(mutinf_res_matrix)
    tenth_kT_mutinfs[tenth_kT_mutinfs < 0.25 * (8.314*300/(4.184*1000)) ] = 0

    output_matrix(resfile_fn+"_tenth_kT_mutinfs.txt",tenth_kT_mutinfs,name_num_list,name_num_list)
    output_matrix(resfile_fn+"_quarter_kT_mutinfs.txt",quarter_kT_mutinfs,name_num_list,name_num_list)
    output_matrix(resfile_fn+"_half_kT_mutinfs.txt",half_kT_mutinfs,name_num_list,name_num_list)
    #now, grab distances for values that pass

    half_kT_CA_dist_matrix = array(CA_dist_matrix)
    half_kT_CA_dist_matrix[half_kT_mutinfs == 0] = 0
    output_matrix(resfile_fn+"_half_kT_CA_dist_matrix.txt",half_kT_CA_dist_matrix,name_num_list,name_num_list)

    quarter_kT_CA_dist_matrix = array(CA_dist_matrix)
    quarter_kT_CA_dist_matrix[quarter_kT_mutinfs == 0] = 0
    output_matrix(resfile_fn+"_quarter_kT_CA_dist_matrix.txt",quarter_kT_CA_dist_matrix,name_num_list,name_num_list)

    tenth_kT_CA_dist_matrix = array(CA_dist_matrix)
    tenth_kT_CA_dist_matrix[tenth_kT_mutinfs == 0] = 0
    output_matrix(resfile_fn+"_tenth_kT_CA_dist_matrix.txt",tenth_kT_CA_dist_matrix,name_num_list,name_num_list)

    filtered_mutinfs = array(mutinf_res_matrix)
    filtered_mutinfs[CA_dist_matrix < options.distance] = 0
    output_matrix(resfile_fn+"_CA_dist_filtered"+str(options.distance)+"_matrix.txt",filtered_mutinfs,name_num_list,name_num_list)

    half_kT_CA_dist_matrix_filtered = array(half_kT_mutinfs)
    half_kT_CA_dist_matrix_filtered[half_kT_mutinfs == 0] = 0
    half_kT_CA_dist_matrix_filtered[CA_dist_matrix < options.distance] = 0
    output_matrix(resfile_fn+"_half_kT_CA_dist_filtered"+str(options.distance)+"_matrix.txt",half_kT_CA_dist_matrix_filtered,name_num_list,name_num_list)

    quarter_kT_CA_dist_matrix_filtered = array(quarter_kT_mutinfs)
    quarter_kT_CA_dist_matrix_filtered[quarter_kT_mutinfs == 0] = 0
    quarter_kT_CA_dist_matrix_filtered[CA_dist_matrix < options.distance] = 0
    output_matrix(resfile_fn+"_quarter_kT_CA_dist_filtered"+str(options.distance)+"_matrix.txt",quarter_kT_CA_dist_matrix_filtered,name_num_list,name_num_list)

    tenth_kT_CA_dist_matrix_filtered = array(tenth_kT_mutinfs)
    tenth_kT_CA_dist_matrix_filtered[tenth_kT_mutinfs == 0] = 0
    tenth_kT_CA_dist_matrix_filtered[CA_dist_matrix < options.distance] = 0
    output_matrix(resfile_fn+"_tenth_kT_CA_dist_filtered"+str(options.distance)+"_matrix.txt",tenth_kT_CA_dist_matrix_filtered,name_num_list,name_num_list)


#END

