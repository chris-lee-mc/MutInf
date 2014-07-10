#!/usr/bin/python

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

#from dihedral_mutent.py import *



def read_res_matrix(myfilename):  #from dihedral_mutent.py 
   rownames = []
   colnames = []
   myfile = open(myfilename,'r')
   inlines = myfile.readlines()
   myfile.close()
   res = inlines[0].split()
   mymatrix = zeros((int(len(inlines[1:])), int(len(res))),float64)
   #print mymatrix.shape
   for myname_num in res:
       colnames.append(myname_num)
   #print colnames
   #print len(colnames)
   for row_num in range(int(len(inlines[1:]))):
       thisline = inlines[row_num + 1]
       thislinedata = thisline.split()
       thisname = thislinedata[0]
       res_num = int(floor(row_num))
       thislinenums = map(float, thislinedata[1:]) #does this need to be float64 or another double precision thingy?
       #print thislinenums
       thislinearray = array(thislinenums,float64)
       #print thislinearray.shape
       rownames.append(thisname)
       for col_num in range(len(colnames)):
           #print "name: "+str(thisname)+" chi: "+str(row_chi)+ " row_num: "+str(row_num)+" row_chi: "+str(row_chi)+ " col_num: "+str(col_num)+" col_chi: "+str(col_chi)+"\n"
           mymatrix[res_num,col_num] = float64(thislinearray[col_num])
   #print rownames
   return mymatrix, rownames, colnames


def run_communities(options, mylambda):
    try:
        fil = open(options.filename,'r')
    except:
        print "cannot open input file "+str(options.filename)
        sys.exit(1) 
    inlines = fil.readlines()
    fil.close()
    
    try:
       logfile = open("log.txt",'w')
    except:
        print "cannot open output file "+str(options.outfile)
        sys.exit(1) 

    ## OPEN COMMUNITIES OUTPUT FILE ##
    try:
        outfile = open(options.outfile+".txt",'w')
    except:
        print "cannot open output file "+str(options.outfile)
        sys.exit(1) 


    

    maxclusternum = -1 #maximum cluster num
    reslist = [] #list of residues


    mycompiler = "gcc"
    
    print "COMMANDS: ", " ".join(sys.argv)
    print "user communities: "
    print options.user_communities 
    #read in data
    mi_matrix, rownames, colnames = read_res_matrix(options.filename)
    mi_matrix += 1e-8
    mi_matrix_old = mi_matrix[:,:] #copy

    contacts, rownames_contacts, colnames_contacts = read_res_matrix(options.matrix)
    
    coords = genfromtxt(options.structure)

    num_res = len(rownames)

    #filter mutual information matrix
    mi_matrix[mi_matrix < options.cutoff] = 0 
    
    #set mutinf = 1 for equal weights
    if options.equal_weights == "yes":
       mi_matrix[:,:] = 1.0

    #make edge list
    print "making edge list"
    logfile.write("making edge list\n")
    edgefile=open(options.edgelist,'w')
    for i in range(num_res):
        for j in range(num_res):
            if (contacts[i,j] > 0 or options.nocontacts == "yes") and mi_matrix[i,j] > 0:
                        edgefile.write("%6i  %6i  %7.4f\n" % (i,j, mylambda(mi_matrix[i,j])  ))
    edgefile.close()
    save_figs = True
    # Load graph from edge list file
    print "reading edge list, making graph"
    logfile.write("reading edge list, making graph \n")
    G = nx.read_weighted_edgelist(options.edgelist, nodetype=int)

    # Perform Girvan Newman algorithm
    print "Girvan Newman clustering"
    logfile.write("Girvan Newman clustering \n")
    if(options.iterations == None):
       itr = int(len(G.edges()) - 9)    # Number of iterations
    else:
       itr = int(options.iterations)
    #itr = int(len(G.edges()) / 10)    # Number of iterations
       
    if options.user_communities==None or options.user_communities=="no":
       max_com_graph, l_com_graph, modularity = gn.Girvan_Newman_algorithm(G, itr,  community_graph=None,
                           target_num_communities=None, logfile=logfile)
       pylab.figure()
       pylab.plot(modularity)
       print 'Maximum modularity = %7.4f' %max(modularity)
       logfile.write('Maximum modularity = %7.4f \n' %max(modularity))
       pylab.xlabel('Iteration')
       pylab.ylabel('Modularity Q-value')
       pylab.title('Modularity')
       if save_figs:
          pylab.savefig('%s_modularity.pdf' % options.prefix, format='pdf')
    
    
    # Show the Original Graph
    pylab.figure()
    ax = pylab.gca()
    positions = nx.layout.spring_layout(G)
    nx.draw(G, ax=ax, pos=positions)
    pylab.title('Network')
    if save_figs:
	pylab.savefig('%s_network.pdf'%options.prefix, format='pdf')

    # Show the community structure graph (all nodes, colored by community)
    pylab.figure()
    if options.user_communities=="no" or  options.user_communities==None:
       communities = gn.find_communities(max_com_graph)
    else:
       (graph, communities, colors) = gn.unpickle_data("communities.pickle")

    print "communities:"
    print communities
    #logfile.write("communities")
    #logfile.write(communities)
    logfile.write("\n")
    gn.draw_communities(G, communities, pos=positions)
    pylab.title('Network with residues colored by Community')
    if save_figs:
	pylab.savefig('%s_community_structure.pdf'%options.prefix, format='pdf')

    # Show coarsegrained community graph 
    # Note: size of nodes proportional to the number of residues
    #       width of edges proportional to the number of shortest paths 
    #       between communities
    pylab.figure()
    cg_graph = gn.coarsegrain_communities(G, communities)

    # Make positions of nodes match up to view in VMD
    # expecting structure file to have to trailing stuff, just ATOM records, one per residue used in mutinf
    
    vmdcolors = [(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1),(0,0,1,1),(1,0,0,1),(.35,0.35,.35,1),(1,.5,0,1),(1,1,0,1),(0.5,0.5,0.2,1),(.45,.0,.9,1),(0,1,0,1),(1,1,1,1),(1,.6,.6,1),(.25,.75,.75,1),(0.65,0,.65,1),(0.5,0.9,0.4,1),(0.9,0.4,0.7,1),(0.5,0.3,0,1),(0.5,0.5,0.75,1) ]
    pos={}
    colors=[]
    for i in range(len(communities)):
        pos.update({i:(mean(coords[communities[i],6]),mean(coords[communities[i],7]))})
        colors.append(vmdcolors[i])

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


    gn.print_coarsegrain(cg_graph, myfilename="cg_edgeweights.txt")
    gn.draw_coarsegrain(cg_graph,pos=pos,colors=colors,node_factor=10.0, edge_factor=50.)
    pylab.title('Community Network Coarsegrain')
    if save_figs:
        pylab.savefig('%s_community_structure_cg.pdf'%options.prefix, format='pdf')
    
    #repalce with edgeweights
    pylab.figure()
    cg_options = cg_edgeweights_to_matrix_for_dendro_default_options()
    cg_options.mutinf = options.filename
    cg_options.contacts = options.matrix
    #grab mutual information between communities
    mutinf_between_communities_contacts_filtered = do_cg_edgeweights_from_mutinf_matrix(cg_options)
    gn.replace_weights_and_draw_coarsegrain(cg_graph, mutinf_between_communities_contacts_filtered, pos=pos,colors=colors,node_factor=10.0, edge_factor=0.005)
    gn.print_coarsegrain(cg_graph, myfilename="cg_mutinf_edgeweights.txt" )
    pylab.title('Community Network Coarsegrain: Mutual Information Edge Weights')
    if save_figs:
        pylab.savefig('%s_community_structure_cg_mutinf_edges.pdf'%options.prefix, format='pdf')



    ## READ COMMUNITIES FILE ##    
    try:
        outfile = open(options.outfile+".txt",'r')
    except:
        print "cannot open output file "+str(options.outfile+".txt")
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
    if(options.display_structure == None):
       options.display_structure = options.structure
    outfile = open(options.outfile+".pml", 'w')
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





    ## show graphics
    #pylab.ion()
    #pylab.show()

    #raw_input("Press ENTER to exit")

    # Save the communities classification to a pickle file
    gn.pickle_data(G, communities)

    #vmdtemp=open('temp_'+options.prefix+'.tcl','w')
    #vmdtemp.write("mol new "+options.structure+" waitfor all\n")

    #for i in range(len(communities)):
    #    val=i+1
    #    #if (val==6):
    #    #    val=26
    #    for j in range(len(communities[i])):
    #                vmdtemp.write("[atomselect top \"residue %i\"] set beta %i\n"%(communities[i][j],val))
    #vmdtemp.write("[atomselect top all] writepdb %s.pdb\nexit\n" % options.prefix)
    #vmdtemp.close()
    #os.system("vmd -dispdev text -e temp_"+options.prefix+"_.tcl")
    # Enter debugger
    # import pdb
    # pdb.set_trace()




if __name__ == "__main__":
    
    usage="%prog [-t traj1:traj2] [-x xvg_basedir] resfile [simulation numbers to use]  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    parser.add_option("-s", "--structure", default=None, type="string", help="pdb file of atoms to use for residue locations")
    parser.add_option("-f", "--filename", default=None, type="string", help="mutual information matrix")
    parser.add_option("-o", "--outfile", default="communities", type="string", help="prefix for output files")
    parser.add_option("-b", "--begin", default=0, type=int, help="first residue offset")
    parser.add_option("-m", "--matrix", default=None, type="string", help="matrix for contacts")
    parser.add_option("-e", "--edgelist", default="edgelist.dat",type="string", help="output edgelist data file")
    parser.add_option("-c", "--cutoff", default=0.5, type="float", help="cutoff value for edges")
    parser.add_option("-p", "--prefix",default="mutinf", type="string", help="prefix for output filenames")
    parser.add_option("-i", "--iterations",default=None, type="int", help="number of iterations")
    parser.add_option("-d", "--display_structure",default=None, type="string", help="full structure for pymol display")
    parser.add_option("-x", "--nocontacts",default="no", type="string", help="whether to use contacts matrix or not ")
    parser.add_option("-U", "--user_communities",default="no",type="string",help="pickle file of previous communities run")
    parser.add_option("-e", "--equal_weights",default="no",type="string", help="whether or not to use mutinf-based edge weights or equal weights")
    
    ## READ INPUT FILE ##    
    (options,args)=parser.parse_args()
    
    print "COMMANDS: ", " ".join(sys.argv)
    mylambda = lambda(x): -.5*log(1-exp(-2*x/3)) #lambda function for converting matrix to edge weights
    run_communities(options, mylambda)
