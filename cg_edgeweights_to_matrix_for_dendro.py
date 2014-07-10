import re, os, sys
from optparse import OptionParser
from numpy import *
import subprocess
from communities import *

class ResCluster :
    members = []  
    selection_text = None #selection text for pymol
    number = 0 #cluster number
    intra_cluster_variance = 0
    extra_cluster_variance = 0
    indexlist = []
    indexlist_complement = []
    def __init__(self,number,firstcluster):
        self.number = number 
        self.members = []
        self.firstcluster = firstcluster
    def __str__(self):
        if(self.selection_text == None):
            self.selection_text = ""
            firstone = self.firstcluster
            for member in self.members:
                if firstone !=1:
                    self.selection_text +=","
                    firstone = 0
                self.selection_text += "(resi "+str(member.number)+" "
                if member.tag == "":
                    self.selection_text += " and (name N,H,CA,C,O) "
                if member.tag == "S":
                    self.selection_text += " and (not name N,H,CA,C,O) "
                if member.chain != None:
                    self.selection_text += " and (chain "+str(member.chain)+" ) "
                self.selection_text += " )"
        return "sele cluster"+str(self.number)+", "+"( "+self.selection_text+")\n color "+str(int(self.number+1))+", cluster"+str(self.number)+"\n"

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

# output the elements of a matrix in string formatting, optionally zeroing the diagonal terms
def output_matrix(myfilename,mymatrix,rownames,colnames, zero_diag=False):
   myfile = open(myfilename,'w')

   for col_num, col_name in zip(range(len(colnames)), colnames):
      myfile.write(str(col_name) + " ")
   myfile.write("\n")
   for row_num, row_name in zip(range(len(rownames)), rownames):
      myfile.write(str(row_name) + " ")
      for col_num, col_name in zip(range(len(colnames)), colnames):
         if col_num == row_num and zero_diag:
            myfile.write(str(0))
         else:
            #print row_num, col_num, mymatrix[row_num,col_num]
            myfile.write(str(mymatrix[row_num,col_num]))
         myfile.write(" ")
      myfile.write("\n")
   myfile.close()


class Res:
    name = "ALA"
    number = 0
    clusternum = -1 #not assigned to a cluster
    tag = "" #default to mainchain
    chain = None #default to no chain
    def __init__(self,name,number,tag,clusternum,chain=None):
        self.name = name
        self.number = number
        self.tag = tag
        self.clusternum = clusternum
        self.chain = chain
    def __str__(self):
        return str(self.name)+" "+str(self.number)+" "+str(self.tag)+" "+str(self.clusternum)+" "+str(self.chain)
    def matrix_header(self):
        if(self.chain == None):
            return str(self.name)+str(self.number)+str(self.tag)
        else:
            return str(self.name)+str(self.number)+str(self.tag)+str(self.chain)
    def __eq__(self, other):
        return ( self.name == other.name and self.number == other.number and self.clusternum == other.clusternum and self.chain == other.chain and self.tag == other.tag ) 

#############################################################################################################################################################################

def do_cg_edgeweights_from_mutinf_matrix(options):

    try:
        Rfile = open(options.Rfile,'w')
    except:
        print "cannot open R command file to be created "+str(options.Rfile)
        sys.exit(1) 
    
    if options.mutinf != None:
        mutinf_matrix, rownames, colnames = read_res_matrix(options.mutinf)

    if options.dist_variance != None:
        dist_variance_matrix, rownames2, colnames2 = read_res_matrix(options.dist_variance)
        print "shape of dist variance matrix: "+str(shape(dist_variance_matrix))

    print "len rownames: "+str(len(rownames))
    print "len colnames: "+str(len(colnames))
    print "mutinf matrix shape: "+str(mutinf_matrix.shape)
    if options.contacts != None:
        try:
            contacts, rownames_contacts, colnames_contacts = read_res_matrix(options.contacts)
        except:
            print "cannot open contacts file: "+str(options.contacts)
        print contacts
    
    #########################################################################


    ## Read input file if a matrix is not provided, or cluster file that was just made by R if a matrix was provided ## 
    

    
    sorted_communities = options.filename
    sorted_communities.replace('.txt','')
    sorted_communities += "_sorted.txt"
    p = subprocess.Popen("cp -p "+str(options.filename)+" "+str(options.filename)+"_bak.txt"+" ; sort -n -k 3 -k 1 "+str(options.filename)+" > "+str(sorted_communities), shell=True)
    os.waitpid(p.pid, 0)
    

    try:
        fil1 = open(options.filename,'r')
    except:
        print "cannot open input file "+str(options.filename)
    try:
        fil2 = open(sorted_communities,'r')
    except:
        print "cannot open input file "+str(sorted_communities)
        sys.exit(1) 
    inlines1 = fil1.readlines()
    fil1.close()
    inlines2 = fil2.readlines()
    fil2.close()

    maxclusternum = -1 #maximum cluster num
    reslist = [] #list of residues
    reslist_sorted = []

    #read communities_sorted.txt file
    

    #myline = inlines[-1]
    #matchline=re.compile(r'([0-9]*)\s*([A-Z][A-Z,0-9][A-Z,0-9])([0-9]+)([S]*)\s+([0-9]+)')
    #matches = matchline.match(myline)
    #if matches.group(2) != None:
    #        name = matches.group(2)
    #        if matches.group(3) != None:
    #            number = int(matches.group(3))
    #            if matches.group(4) == 'S':
    #                tag = matches.group(4)
    #            else:
    #                tag = ''
    #            if matches.group(5) != None:
    #                clusternum = int(matches.group(5))
    #                last_community = clusternum

    
    community_list = []  #
    templist = []
    communities_lengths = []
    for line in inlines1[0:]:
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

    this_community = 0
    for line in inlines2[0:]:
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
                    if clusternum > this_community:
                        communities_lengths.append(len(templist))
                        community_list.append(templist)
                        templist = [] #start building up next community 
                        this_community = clusternum
                        print "community boundary: "+str(line)
                    newres = Res(name,number,tag,clusternum) #offset residue number by option
                    reslist_sorted.append(newres)
                    templist.append(newres)
                    if(clusternum > maxclusternum):
                        maxclusternum = clusternum
    #last iteration
    communities_lengths.append(len(templist))
    community_list.append(templist)
    this_community = clusternum
    print "community boundary: "+str(line)
                        
                        
    # here's where we build up the new mutual information matrix 
    print "max community number: "+str(maxclusternum)
    print "length of reslist: "+str(len(reslist))
    print "communities lengths:"+str(communities_lengths)
    communities_sorted_mutinf = mutinf_matrix.copy()
    rownames_new = []
    colnames_new = []
    mutinf_contacts_filtered = zeros((len(reslist), len(reslist)), float64)
    mutinf_between_communities = zeros((maxclusternum + 1, maxclusternum + 1), float64)
    mutinf_between_communities_contacts_filtered = zeros((maxclusternum + 1, maxclusternum + 1), float64)
    communities_sorted_mutinf_blocks = zeros((len(reslist), len(reslist)), float64)
    communities_sorted_mutinf_contacts_filtered = zeros((len(reslist), len(reslist)), float64)
    communities_sorted_mutinf_contacts_filtered_blocks = zeros((len(reslist), len(reslist)), float64)
    avg_distance_variances_within_and_between_communities = zeros((maxclusternum + 1, maxclusternum + 1), float64)
    
    
    for i, myres1 in zip(range(len(reslist_sorted)),reslist_sorted):
        for ik, ikres1 in zip(range(len(reslist)), reslist):
            if myres1 == ikres1: #another way of getting index 
                myindex1 = ik
        rownames_new.append(rownames[myindex1] )
        colnames_new.append(rownames[myindex1] )
        for j, myres2 in zip(range(len(reslist_sorted)),reslist_sorted):
            myindex2 = reslist.index(myres2)
            for jk, jkres2 in zip(range(len(reslist)), reslist):
                if myres2 == jkres2:
                    myindex2 = jk
            communities_sorted_mutinf[i,j] = mutinf_matrix[myindex1,myindex2] #reorders matrix
            if options.contacts != None:
                communities_sorted_mutinf_contacts_filtered[i,j] = contacts[myindex1,myindex2] *  mutinf_matrix[myindex1,myindex2] 
                mutinf_between_communities_contacts_filtered[myres1.clusternum, myres2.clusternum] +=  contacts[myindex1,myindex2] * mutinf_matrix[myindex1,myindex2] #/ (1.0 * len(community_list[i]) * len(community_list[j]))
                mutinf_contacts_filtered[i,j] = contacts[i,j] * mutinf_matrix[i,j]
            if options.dist_variance != None:
                #print "dist matrix element: "+str( dist_variance_matrix[myindex1,myindex2])
                #print "clusternum1 :"+str(myres1.clusternum)
                #print "clusternum2 :"+str(myres2.clusternum)
                                                  
                avg_distance_variances_within_and_between_communities[myres1.clusternum, myres2.clusternum] += dist_variance_matrix[myindex1,myindex2] / (1.0 * communities_lengths[myres1.clusternum] * communities_lengths[myres2.clusternum])
            mutinf_between_communities[myres1.clusternum, myres2.clusternum] +=  mutinf_matrix[myindex1,myindex2] / (1.0 * communities_lengths[myres1.clusternum] * communities_lengths[myres2.clusternum]) #average mutinf btw communities
             
    for i, myres1 in zip(range(len(reslist_sorted)),reslist_sorted):
        myindex1 = reslist.index(myres1) 
        for j, myres2 in zip(range(len(reslist_sorted)),reslist_sorted):
            communities_sorted_mutinf_blocks[i,j] = mutinf_between_communities[myres1.clusternum, myres2.clusternum]
            if options.contacts != None:
                communities_sorted_mutinf_contacts_filtered_blocks[i,j] = mutinf_between_communities_contacts_filtered[myres1.clusternum, myres2.clusternum] 
                
                
    output_matrix(options.prefix+options.outfile,communities_sorted_mutinf,rownames_new,colnames_new)
    

    rownames_communities = []
    for i in range(maxclusternum + 1):
        rownames_communities.append("C"+str(i+1))
    colnames_communities = rownames_communities
    output_matrix(str(options.prefix)+options.communities_mutinf,  mutinf_between_communities, rownames_communities, colnames_communities)
    output_matrix(str(options.prefix)+options.communities_mutinf+"_avg_distance_variances.txt", avg_distance_variances_within_and_between_communities , rownames_communities, colnames_communities)

    output_matrix(str(options.prefix)+options.communities_mutinf+"_blocks.txt",  communities_sorted_mutinf_blocks, rownames_new, colnames_new)
    output_matrix(str(options.prefix)+"mutinf_unclustered_contacts_filtered.txt", mutinf_contacts_filtered , rownames, colnames)
    ### NEED TO ADD ABOVE TO R SCRIPT, and HOOKS TO SPECIFY DISTANCE VARIANCE MATRICES, this will tell us if communities are semi-rigid elements ###
    output_matrix(str(options.prefix)+options.communities_mutinf+"_contacts_filtered.txt",  communities_sorted_mutinf_contacts_filtered, rownames_new, colnames_new)
    output_matrix(str(options.prefix)+options.communities_mutinf+"_contacts_filtered_blocks.txt",  communities_sorted_mutinf_contacts_filtered_blocks, rownames_new, colnames_new)

    num_residues = len(reslist)
    
    R_code1 = """
       library(marray)
       library(fields)  
       library(cluster)
       """
    Rfile.write(R_code1)
    #Rfile.write("prec <- read.table(\""+str(options.matrix)+"\")\n")
    #Rfile.write("precdist <- read.table(\""+str(distfile)+"\")\n")
    Rfile.write("mutinf_unclustered <- read.table(\""+str(options.mutinf)+"\")\n")
    Rfile.write("new_mutinf <- read.table(\""+str(options.prefix)+str(options.outfile)+"\")\n")
    Rfile.write("mutinf_unclustered_contacts_filtered <- read.table(\""+str(options.prefix)+"mutinf_unclustered_contacts_filtered.txt"+"\")\n")
    Rfile.write("communities_mutinf <- read.table(\""+str(options.prefix)+str(options.communities_mutinf)+"\")\n")
    Rfile.write("communities_mutinf_blocks <- read.table(\""+str(options.prefix)+str(options.communities_mutinf+"_blocks.txt")+"\")\n")
    Rfile.write("communities_mutinf_contacts_filtered <- read.table(\""+str(options.prefix)+str(options.communities_mutinf+"_contacts_filtered.txt")+"\")\n")
    Rfile.write("communities_mutinf_contacts_filtered_blocks <- read.table(\""+str(options.prefix)+str(options.communities_mutinf+"_contacts_filtered_blocks.txt")+"\")\n")
    Rfile.write("pdf(\""+str(options.prefix)+"mutinf_communities_reordered.pdf""\", width=96,height=96)\n")
    mydict = {'mutinf_unclustered': options.prefix+"mutinf_unclustered.pdf" ,'mutinf_unclustered_contacts_filtered' : str(options.prefix)+"mutinf_unclustered_contacts_filtered.pdf", 'mutinf_between_communities': options.prefix+"mutinf_between_communities.pdf", 'mutinf_between_communities_contacts_filtered_blocks': options.prefix+"mutinf_between_communities_contacts_filtered_blocks.pdf", 'mutinf_between_communities_blocks': options.prefix+"mutinf_between_communities_blocks.pdf",'mutinf_between_communities_contacts_filtered': options.prefix+"mutinf_between_communities_contacts_filtered.pdf", 'colorbar_mutinf': options.prefix+"colorbar_mutinf.pdf", 'colorbar_mutinf_between_communities': options.prefix+"colorbar_mutinf_between_communities.pdf"}
    R_code2 = """
               library(marray)

               #heatmap(as.matrix(new_mutinf), col=maPalette(low="white",mid="blue",high="red", k=50), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)
               heatmap(as.matrix(new_mutinf), col=maPalette(low="white",mid="gray",high="red", k=50), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)
               pdf("%(mutinf_unclustered)s", width=96, height=96)
               heatmap(as.matrix(mutinf_unclustered), col=maPalette(low="white",mid="gray",high="red", k=50), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)
               pdf("%(mutinf_unclustered_contacts_filtered)s", width=96, height=96)
               heatmap(as.matrix(mutinf_unclustered_contacts_filtered), col=maPalette(low="white",mid="gray",high="red", k=50), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)              
  
               #image(as.matrix(new_mutinf),col = maPalette(low = "white", high = "black", k =50), add = FALSE, xaxt= "n", yaxt = "n",xlab="Residue i",ylab="Residue j")
               #axis(2,labels=c(0,20,40,60,80,100,120,132),at=c(1-0,1-20/132,1-40/132,1-60/132,1-80/132,1-100/132,1-120/132,1-132/132))
               #axis(1,labels=c(0,20,40,60,80,100,120,132),at=c(0,20/132,40/132,60/132,80/132,100/132,120/132,132/132))
               #rect(72/132,1-46/132,86/132,1-29/132,border="red")
               #rect(29/132,1-46/132,46/132,1-29/132,border="red")
               #rect(72/132,1-86/132,86/132,1-72/132,border="red")
               #rect(72/132,1-6/132,86/132,1-4/132,border="blue")

 dev.off()

               pdf("%(mutinf_between_communities)s")
               #heatmap(as.matrix(communities_mutinf), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, oldstyle=FALSE,symm = TRUE)
               heatmap(as.matrix(communities_mutinf), col=maPalette(low="white",mid="gray",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, oldstyle=FALSE,symm = TRUE)
               pdf("%(mutinf_between_communities_contacts_filtered_blocks)s",width=96,height=96)
               #heatmap(as.matrix(communities_mutinf_contacts_filtered_blocks), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, oldstyle=FALSE,symm = TRUE)
               heatmap(as.matrix(communities_mutinf_contacts_filtered_blocks), col=maPalette(low="white",mid="gray",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA, oldstyle=FALSE,symm = TRUE)
               pdf("%(mutinf_between_communities_blocks)s",  width=96,height=96)
               #heatmap(as.matrix(communities_mutinf_blocks), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA,  cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)
               heatmap(as.matrix(communities_mutinf_blocks), col=maPalette(low="white",mid="gray",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA,  cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)

               pdf("%(mutinf_between_communities_contacts_filtered)s",  width=96,height=96)
               #heatmap(as.matrix(communities_mutinf_contacts_filtered), col=maPalette(low="white",mid="blue",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA,  cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)
               heatmap(as.matrix(communities_mutinf_contacts_filtered), col=maPalette(low="white",mid="gray",high="red"), add = FALSE, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", Rowv = NA, Colv = NA,  cexRow=0.8, cexCol=0.8, oldstyle=FALSE,symm = TRUE)

               pdf("%(colorbar_mutinf)s",width=6,height=2)
               #maColorBar(seq(min(new_mutinf),max(new_mutinf),by=(max(new_mutinf)-min(new_mutinf))/50), horizontal=TRUE, col = maPalette(low = "white", mid="blue", high = "red", k =50))
               maColorBar(seq(min(new_mutinf),max(new_mutinf),by=(max(new_mutinf)-min(new_mutinf))/50), horizontal=TRUE, col = maPalette(low = "white", mid = "gray", high = "red", k =50))

               pdf("%(colorbar_mutinf_between_communities)s",width=6,height=2)
               #maColorBar(seq(min(communities_mutinf),max(communities_mutinf),by=(max(communities_mutinf)-min(communities_mutinf))/50), horizontal=TRUE, col = maPalette(low = "white", mid="blue", high = "red", k =50))
               maColorBar(seq(min(communities_mutinf),max(communities_mutinf),by=(max(communities_mutinf)-min(communities_mutinf))/50), horizontal=TRUE, col = maPalette(low = "white", mid="gray", high="red", k =50))
               
#heatmap(x=as.matrix((S1M47_mutent)),zlim=c(0,max(S1M47_mutent)*1.0), col = maPalette(low = "white", mid = "blue", high = "red", k =50),     add = FALSE, xaxs = "i", yaxs = "i",xaxt= "n", yaxt = "n", reorderfun = #function(d,w) rev(reorder(d,w)),revC=TRUE, cexRow=0.8,  cexCol=0.8,      oldstyle = FALSE,symm=TRUE )
#dev.off()
          """ % mydict
    
    Rfile.write(R_code2)
    #Rfile.write("write.table(fit$cluster, \""+str(options.filename)+"\", quote=FALSE, sep='\t')")
    Rfile.write("\n")
    Rfile.write("\n")
    Rfile.close()
    
    ### Run R to calculate clusters using multidimensional scaling  #######
    
    print "Running R on "+str(options.Rfile)
    p = subprocess.Popen("cat "+str(options.Rfile)+" | R --no-save", shell=True)
    os.waitpid(p.pid, 0)
    
    return mutinf_between_communities
    #return mutinf between communities filtered by contacts for further analysis, if desired
    #return mutinf_between_communities_contacts_filtered

def cg_edgeweights_to_matrix_for_dendro_default_options():

    class run_options:
        filename = "communities.txt"
        outfile = "mutinf_communities_reordered.txt"
        begin = 0
        Rfile = "heatmap_communities.txt"
        mutinf = None
        communities_mutinf = "mutinf_bewteen_communities.txt"
        mutinf_output = "mutinf_unclustered.txt"
        contacts = None
        prefix = ""
        dist_variance=None

    my_options = run_options()
    return my_options

if __name__ == "__main__":
    usage="%prog -f communities.txt -t mypdb_mutinf_res_sum_0diag.txt "
    parser=OptionParser(usage)
    ### options for k-means and multidimensional scaling
    #parser.add_option("-s", "--structure", default=None, type="string", help="pdb file")
    parser.add_option("-f", "--filename", default="communities.txt", type="string", help="space-delimeted text file with three columns: number  residue name/number/tag  cluster_number")
    parser.add_option("-o", "--outfile", default="mutinf_communities_reordered.txt", type="string", help="filename for output mutual information matrix ordered by communities")
    parser.add_option("-b", "--begin", default=0, type=int, help="first residue offset")
    parser.add_option("-r", "--Rfile", default="heatmap_communities.txt",type="string", help="R commands file to be created")
    parser.add_option("-t", "--mutinf", default=None, type="string", help="mutual information matrix filename")
    parser.add_option("-m", "--communities_mutinf", default="mutinf_bewteen_communities.txt", type="string", help="mutual information between communities matrix filename")
    parser.add_option("-c", "--contacts", default=None, type="string", help="matrix for contacts") 
    parser.add_option("-p", "--prefix", default="",type="string", help="prefix for output")
    parser.add_option("-v", "--dist_variance", default=None, type="string", help="filename for distance variances")
    ## Generate cluster list using multidimensional scaling and kmeans in R ##
    (options,args)=parser.parse_args()
    print "options"
    print options
    #distfile = options.matrix
    #distfile = distfile.replace('dist_variance', 'dist')
    
    mutinf_between_communities = do_cg_edgeweights_from_mutinf_matrix(options)
