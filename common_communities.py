#

from numpy import *
import os
import sys
import cPickle
import girvan_newman as gn
import networkx as nx
import pylab
from optparse import OptionParser
from cluster_list_to_pml import *
from subprocess import call
from cg_edgeweights_to_matrix_for_dendro import *


def list_to_dict(li): 
    dct = {} 
    dct = {key: value for (key, value) in zip(range(len(li)),li) }
    return dct 


vmdcolors = [(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1) ]

save_figs = True

#3FJQ_twoMg

def renumber_communities(options):
    
    community_maps = []
    community_data = []
    minlen = 99
    i = 0
    cutoffs = [8.0, 9.0, 10.0]
    which_minlen = zeros((4),int8)
    print "which_minlen shape: "+str(shape(which_minlen))
    k = 0
    mutinf_filenames = ["3FJQ_twoMg/mutinf_3FJQ_anton_carts_hotatoms_1.0/3FJQ_cart_hotatoms.reslist-nsims6-structs7366-bin30_bootstrap_avg_mutinf_res_sum_0diag.txt",  "3FJQ_atp_Mg2_HIP87/mutinf_3FJQ_anton_carts_hotatoms_1.0/3FJQ_cart_Mg2_hotatoms.reslist-nsims6-structs7936-bin30_bootstrap_avg_mutinf_res_sum_0diag.txt", "1CMK_atp_Mg2_mouse_HIP87/mutinf_1CMK_atp_Mg2_carts_hotatoms_1.0/1CMK_cart_Mg2_hotatoms.reslist-nsims6-structs10834-bin30_bootstrap_avg_mutinf_res_sum_0diag.txt", "1CMK_apo/mutinf_1CMK_apo_HIP87_carts_hotatoms_1.0/1CMK_apo_cart_hotatoms.reslist-nsims6-structs3311-bin30_bootstrap_avg_mutinf_res_sum_0diag.txt"]

    structure_list = ["3FJQ_atp_Mg2_HIP87/hotatoms1g.pdb", "3FJQ_atp_Mg2_HIP87/hotatoms1g.pdb", "3FJQ_atp_Mg2_HIP87/hotatoms1g.pdb" ,"3FJQ_atp_Mg2_HIP87/hotatoms1g.pdb"] # need hotatoms structure file for community centers
    #structure_list = ["3FJQ_atp_Mg2_HIP87/3FJQ_atp_oneMg2_noPKI_maestro_HIP87onlyHIE142_CYM199_anton_amber_nowat.pdb","3FJQ_atp_Mg2_HIP87/3FJQ_atp_oneMg2_noPKI_maestro_HIP87onlyHIE142_CYM199_anton_amber_nowat.pdb","1CMK_atp_Mg2_mouse_HIP87/1CMK_atp_sapiens_HIP87_Mg2_anton_nowat_amber.pdb","1CMK_apo/1CMK_apo_HIP87_anton_nowat_amber.pdb"]

    
    directory_list =  ["3FJQ_twoMg/mutinf_3FJQ_anton_carts_hotatoms", "3FJQ_atp_Mg2_HIP87/mutinf_3FJQ_anton_carts_hotatoms", "1CMK_atp_Mg2_mouse_HIP87/mutinf_1CMK_atp_Mg2_carts_hotatoms",  "1CMK_apo/mutinf_1CMK_apo_HIP87_carts_hotatoms" ]
    
    upper_directory_list = ["3FJQ_twoMg", "3FJQ_atp_Mg2_HIP87", "1CMK_atp_Mg2_mouse_HIP87",  "1CMK_apo"] 

    k = 0
    cutoffs = [8.0, 9.0, 10.0]
    for mydir_prefix in directory_list:
        i = 0
        community_data = []
        for cutoff in [8.0, 9.0, 10.0]:
            
            (graph, communities, colors) = gn.unpickle_data(str(mydir_prefix)+"_"+str(cutoff)+"_1.0/communities.pickle")
            community_data.append(tuple([graph, communities, colors])) #has to be a tuple
            mylen = len(communities)
            print "mylen: "+str(mylen)
            if mylen <= minlen:
                print "k: "+str(k)+" i: "+str(i)
                which_minlen[k] = i
                minlen = mylen
            i += 1
        print "community_data shape "
        print shape(community_data)
        print str(mydir_prefix)+"which minlen: "+str(which_minlen[k])+" cutoff: "+str(cutoffs[which_minlen[k]])
        #print community_data[which_minlen[k]]
        community_maps.append((community_data[which_minlen[k]])[1]) # second element gives the communities
        k += 1
        
    print "community maps shape"
    print shape(community_maps)
    print "community maps shape first item"
    print shape(community_maps[0])
    print "community maps"
    print community_maps

    #map to common communities set 
    (community_maps_new_sets, basis_set) = gn.map_communities(community_maps)
    #
    #print "basis set"
    #print basis_set
    #print "community_maps_new_sets"
    #print community_maps_new_sets
    
    #rebuild as a list
    #print "shape of community_maps_new_sets"
    #print community_maps_new_sets.

    community_maps_new = [] # list(community_maps_new_sets)
    for community_set in community_maps_new_sets:
        community_map = [] 
        for mykey, value in enumerate(community_set):
            print "mykey : "+str(mykey)
            print "value: "
            print value
            print 
            community_map.append(list(community_set[value]))
        community_maps_new.append(community_map)
    #community_maps_new = community_maps   

    print "new community maps first element"
    print community_maps_new[0]
    k = 0
    for mydir_prefix in directory_list:
       mi_matrix, rownames, colnames = read_res_matrix(mutinf_filenames[k])
       G =  nx.read_weighted_edgelist(mydir_prefix+"_"+str(cutoffs[which_minlen[k]])+"_1.0/"+str(options.edgelist), nodetype=int)
       communities = community_maps_new[k]
       # Show coarsegrained community graph 
       # Note: size of nodes proportional to the number of residues
       #       width of edges proportional to the number of shortest paths 
       #       between communities
       pylab.figure()
       cg_graph = gn.coarsegrain_communities(G, communities)

       # Make positions of nodes match up to view in VMD
       # expecting structure file to have to trailing stuff, just ATOM records, one per residue used in mutinf
    
       ## OPEN COMMUNITIES OUTPUT FILE ##
       options.outfile += "_renumbered"   #this is important, to not overwrite originals
       #myprefix = upper_directory_list[k]+"/"+str(cutoffs[which_minlen[k]])+"_1.0/"
       myprefix = upper_directory_list[k]+"/"
       try:
           outfile = open(str(myprefix)+options.outfile+".txt",'w')
       except:
           print "cannot open output file "+str(str(myprefix)+options.outfile+".txt")
           sys.exit(1) 
       
       options.structure = structure_list[k]
       coords = genfromtxt(options.structure)
       
       pos={}
       colors=[]
       for i in range(len(communities)):
           pos.update({i:(mean(coords[communities[i],6]),mean(coords[communities[i],7]))})
           colors.append(vmdcolors[i])
       print pos

       ## Output communities to file
       community_mappings=zeros((len(rownames)),int16)
       for i in range(len(communities)):
           for j in range(len(communities[i])):
               #print "communities:"
               #print communities[i]
               #print "communities i j:"
               #print communities[i][j]
               if(communities[i][j] >= 0):
                   community_mappings[communities[i][j]] = i
       
       
       for i in range(len(rownames)):
           outfile.write(str(i)+"   "+str(rownames[i])+"   "+str(community_mappings[i])+"\n")
       outfile.close()
       
       
       gn.print_coarsegrain(cg_graph, myfilename="%s_cg_edgeweights_renumbered.txt"%options.prefix)
       gn.draw_coarsegrain(cg_graph,pos=pos,colors=colors,node_factor=10.0, edge_factor=50.)
       pylab.title('Community Network Coarsegrain')
       if save_figs:
           pylab.savefig(myprefix+'%s_community_structure_cg_renumbered.pdf'%options.prefix, format='pdf')
       
       #repalce with edgeweights
       pylab.figure()
       cg_options = cg_edgeweights_to_matrix_for_dendro_default_options()
       cg_options.mutinf = mutinf_filenames[k]
       cg_options.filename = str(myprefix)+options.outfile+".txt"
       cg_options.contacts = None
       #grab mutual information between communities
       #mutinf_between_communities = do_cg_edgeweights_from_mutinf_matrix(cg_options)
       #gn.replace_weights_and_draw_coarsegrain(cg_graph, mutinf_between_communities, pos=pos,colors=colors,node_factor=10.0, edge_factor=0.01)
       #gn.print_coarsegrain(cg_graph, myfilename="cg_mutinf_edgeweights.txt" )
       #pylab.title('Community Network Coarsegrain: Mutual Information Edge Weights')
       #if save_figs:
       #    pylab.savefig(myprefix+'%s_community_structure_cg_renumbered_mutinf_edges.pdf'%options.prefix, format='pdf')
       
       
       ## READ COMMUNITIES FILE ##    
       try:
           outfile = open(str(myprefix)+options.outfile+".txt",'r')
       except:
           print "cannot open output file "+str(str(myprefix)+options.outfile+".txt")
           sys.exit(1) 
       inlines = outfile.readlines()
       outfile.close()
       
       maxclusternum = -1 #maximum cluster num
       reslist = [] #list of residues
       
       #read file
       for line in inlines[1:]:
           matchline=re.compile(r'([0-9]*)\s*([A-Z][A-Z,0-9][A-Z,0-9])([0-9]+)([S]*)\s+([0-9]+)')
           #print line
           matches = matchline.match(line)
           #print [matches.group(i) for i in range(5)]
           if matches.group(2) != None:
               name = matches.group(2)
               if matches.group(3) != None:
                   number = int(matches.group(3))
                   if matches.group(4) == 'S':
                       tag = matches.group(4)
                   else:
                       tag = ''
                   if matches.group(5) != None:
                       clusternum = int(matches.group(5))
                       newres = Res(name,number,tag,clusternum) #offset residue number by option
                       reslist.append(newres)
                       if(clusternum > maxclusternum):
                           maxclusternum = clusternum
       
       #write out clusters to pymol file 
       #print "Writing Pymol Session File"
       
       options.display_structure = options.structure
       outfile = open(str(myprefix)+options.outfile+".pml", 'w')
       outfile.write("from pymol import cmd\n")
       outfile.write("cmd.bg_color('white')\n")
       outfile.write("load "+str(options.display_structure)+", system\n")
       clusterlist = []
       #loop over clusters
       for clusternum in range(0,maxclusternum+1):
           thiscluster = []
           thiscluster = ResCluster(clusternum, 0)
           for residue in reslist:
               #print "residue number: "+str(residue.number)+" cluster: "+str(residue.clusternum)
               if residue.clusternum == clusternum:
                   thiscluster.members.append(residue)
           #print "cluster: "+str(clusternum)
           #for mymember in thiscluster.members:
           #    print str(mymember)
           clusterlist.append(thiscluster)
           outfile.write(str(thiscluster)) #output pymol selection line
       
       
       #finish session file
       outfile.write("cmd.show('cartoon'   ,'system')\n")
       #if options.matrix != None:
       #    outfile.write("preset.b_factor_putty('system')"+"\n")  #,_self=cmd"+"\n")
       #    outfile.write("cmd.spectrum('b',selection=('all'),quiet=0")
       outfile.write("cmd.show('sticks','((byres (system))&(!(n;c,o,h|(n. n&!r. pro))))')\n")
       outfile.write("cmd.hide('(all and hydro)')\n")
       outfile.write("alter system, resi=str(int(resi)+"+str(options.begin)+")")
       outfile.close()
       k += 1
       

    
if __name__ == "__main__":
    
    usage="%prog [-t traj1:traj2] [-x xvg_basedir] resfile [simulation numbers to use]  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    parser.add_option("-s", "--structure", default=None, type="string", help="pdb file of atoms to use for residue locations")
    parser.add_option("-o", "--outfile", default="communities", type="string", help="prefix for output files")
    parser.add_option("-b", "--begin", default=0, type=int, help="first residue offset")
    parser.add_option("-e", "--edgelist", default="edgelist.dat",type="string", help="output edgelist data file")
    parser.add_option("-p", "--prefix",default="mutinf", type="string", help="prefix for output filenames")
    
    
    

    (options,args)=parser.parse_args()
    
    print "COMMANDS: ", " ".join(sys.argv)
    renumber_communities(options)

#END
