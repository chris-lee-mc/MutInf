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
        return "sele cluster"+str(self.number+1)+", "+"( "+self.selection_text+")\n color "+str(int(self.number+1))+", cluster"+str(self.number+1)+"\n"

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

if __name__ == "__main__":
    usage="%prog -s mypdb.pdb -f clusters.txt -o clusters.pml"
    parser=OptionParser(usage)
    ### options for k-means and multidimensional scaling
    parser.add_option("-s", "--structure", default=None, type="string", help="pdb file")
    parser.add_option("-f", "--filename", default=None, type="string", help="space-delimeted text file with three columns: number  residue name/number/tag  cluster_number")
    parser.add_option("-o", "--outfile", default="clusters.pml", type="string", help="filename for pymol session file")
    parser.add_option("-b", "--begin", default=0, type=int, help="first residue offset")
    parser.add_option("-m", "--matrix", default=None, type="string", help="matrix for intra/extra cluster average variances")
    parser.add_option("-r", "--Rfile", default="cluster_kmeans_k12.txt",type="string", help="R commands file to be created")
    ### options for community analysis
    parser.add_option("-c", "--contacts", default=None, type="string", help="matrix for contacts") 
    parser.add_option("-u", "--cutoff",default=0.5, type="float", help="cutoff value for edges") 
    parser.add_option("-p", "--prefix",default="mutinf", type="string", help="prefix for output filenames")
    parser.add_option("-i", "--iterations",default=None, type="int", help="number of iterations")
    parser.add_option("-a", "--hotatoms", default="../hotatoms1g.pdb", type="string", help="filename of CA and sidechain locations to use for community analysis")
    parser.add_option("-t", "--mutinf", default=None, type="string", help="mutual information matrix filename")
    parser.add_option("-n", "--noR", default="yes", type="string", help="whether to run R or not (yes/no)")
    parser.add_option("-x", "--nocontacts", default="no", type="string", help="whether to use contacts for filtering mutinf matrix for community analysis")
    parser.add_option("-e", "--equal_weights",default="no",type="string", help="whether or not to use mutinf-based edge weights or equal weights")
    parser.add_option("-U", "--user_communities",default=None,type="string", help="filename for communities pickle")
    #parser.add_option("-d", "--display_structure",default=None, type="string", help="full structure for pymol display")
    ## Generate cluster list using multidimensional scaling and kmeans in R ##
    (options,args)=parser.parse_args()
    print "options"
    print options
    distfile = options.matrix
    distfile = distfile.replace('dist_variance', 'dist')
    
    try:
        Rfile = open(options.Rfile,'w')
    except:
        print "cannot open R command file to be created "+str(options.Rfile)
        sys.exit(1) 
    
    if options.matrix != None:

        variance_matrix, rownames, colnames = read_res_matrix(options.matrix)

    R_code1 = """
       library(marray)
       library(fields)  
       library(cluster)
       """
    Rfile.write(R_code1)
    Rfile.write("prec <- read.table(\""+str(options.matrix)+"\")\n")
    Rfile.write("precdist <- read.table(\""+str(distfile)+"\")\n")

    R_code2 = """
       pdf("distance_variance_vs_distance.pdf",height=12,width=12)
       plot(c(as.matrix(precdist)),c(as.matrix(prec)),xlab="Distance",ylab="Distance Variance")
       dev.off()

       pdf("distance_variance_vs_distance_boxplot.pdf",height=12,width=12)
       bplot.xy(x=c(as.matrix(precdist)),y=c(as.matrix(prec)+diag(nrow(prec))),N=20,xlab="Distance",ylab="Distance Variance")
       dev.off()
       nrow(prec)
       ncol(prec)
       blah = hclust(as.dist(prec),method="ward")
       blah2 = cutree(blah, k=16)
       blah3=as.table(blah2)
       blah4=data.frame(blah3)
       blah4
       # write.table(blah4, "output_clusters.txt", quote=FALSE, sep='\t')

       d <- as.dist(prec) 
       fit <- cmdscale(d,eig=TRUE,k=(nrow(prec)-1)/2) # k is the number of dim
       # fit # view results

       newdata = fit$points
       # plot solution 
       x <- fit$points[,1]
       y <- fit$points[,2]

       plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
         main="Metric	MDS",	 type="n")
       text(x, y, labels = row.names(newdata), cex=.1)

       plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
         main="Metric	MDS",	 type="n")
       text(x, y, labels = row.names(newdata), cex=.1)


       # Determine number of clusters
       wss <- (nrow(prec)-1)*sum(apply(newdata,2,var))
       for (i in 2:14)  {
           ktest <- kmeans(newdata, centers=i, iter=500)
           mysum <- 0
          mysum2 <- 0
          num1 <- 0
          num2 <- 0
          myratio <- 0
          for(j in 1:i)     {
          if(length(which(ktest$cluster == j)) > 1) {
             mysum <- mysum + mean(prec[which(ktest$cluster == j),which(ktest$cluster == j)])
             num1 <- length(which(ktest$cluster == j))
             mysum2 <- mysum2 + mean(prec[which(ktest$cluster == j),which(ktest$cluster != j)])
             myratio <-myratio + mean(prec[which(ktest$cluster == j),which(ktest$cluster == j)]) / mean(prec[which(ktest$cluster == j),which(ktest$cluster != j)])
             num2 <- length(which(ktest$cluster != j))
             }
          }
       #    wss[i] <- (mysum / num1) / (mysum2 / num2)
       #      wss[i] <- (mysum / mysum2)
            wss[i] <- (myratio / num1)
       }

       #wss[i] <- sum(kmeans(newdata,   	 centers=i)$withinss)

       pdf("sum_intra_cluster_distance_variance_vs_clusters.pdf")
       plot(5:14, wss[5:14], type="b", xlab="Number of Clusters",
         ylab="Within groups average intra/extra cluster distance variance ratio")
       dev.off()

       plot(5:14, wss[5:14], type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares")

       # K-Means Cluster Analysis
       fit <- kmeans(newdata, 12, iter=1000) # selected cluster solution

       # get cluster means 
       aggregate(newdata,by=list(fit$cluster),FUN=mean)

       # append cluster assignment
       mydata <- data.frame(newdata, fit$cluster)
       fit$size
       """

    Rfile.write(R_code2)
    Rfile.write("write.table(fit$cluster, \""+str(options.filename)+"\", quote=FALSE, sep='\t')")
    Rfile.write("\n")
    Rfile.write("\n")
    Rfile.close()
    
    ### Run R to calculate clusters using multidimensional scaling  #######
    if options.noR == "yes":
        print "Running R on "+str(options.Rfile)
        p = subprocess.Popen("cat "+str(options.Rfile)+" | R --no-save", shell=True)
        os.waitpid(p.pid, 0)
    
    #########################################################################


    ## Read input file if a matrix is not provided, or cluster file that was just made by R if a matrix was provided ##    
    (options,args)=parser.parse_args()
    try:
        fil = open(options.filename,'r')
    except:
        print "cannot open input file "+str(options.filename)
        sys.exit(1) 
    inlines = fil.readlines()
    fil.close()

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
    print "Writing Pymol Session File "+str(options.outfile)
    outfile = open(options.outfile, 'w')
    outfile.write("from pymol import cmd\n")
    outfile.write("cmd.bg_color('white')\n")
    outfile.write("load "+str(options.structure)+", system\n")
    ## need to fix this if it is to be used ## outfile.write("alter all, resi=str(int(resi -"+str(options.begin)-1")) \n")
    clusterlist = []
    #loop over clusters
    for clusternum in range(1,maxclusternum+1):
        thiscluster = []
        thiscluster = ResCluster(clusternum, 1)
        for residue in reslist:
            #print "residue number: "+str(residue.number)+" cluster: "+str(residue.clusternum)
            if residue.clusternum == clusternum:
                thiscluster.members.append(residue)
        #print "cluster: "+str(clusternum)
        #for mymember in thiscluster.members:
        #    print str(mymember)
        clusterlist.append(thiscluster)
        outfile.write(str(thiscluster)) #output pymol selection line
        
    #if intra-cluster variances over extra-cluster variance averages are to be put on the structure as b-factors
        
        
        #print rownames
        clusternum = 1
        
        for thiscluster in clusterlist:
            clusternum = thiscluster.number
            thiscluster.indexlist = []
            thiscluster.indexlist_complement = range(len(rownames))
            for mymember in thiscluster.members:
                myheader = mymember.matrix_header()
                myindex = rownames.index(myheader) # look for header in rownames
                #print "match "+str(myindex)
                thiscluster.indexlist.append(myindex)
                thiscluster.indexlist_complement.remove(myindex)
            #print "cluster "+str(clusternum)+" members: "+str(len(thiscluster.members))
            #print "within cluster: "+str(len(thiscluster.indexlist))
            #print thiscluster.indexlist
            #for myindex in thiscluster.indexlist:
            #    print rownames[myindex]
            #print "outside cluster: "+str(len(thiscluster.indexlist_complement))
            #print thiscluster.indexlist_complement
            #for myindex in thiscluster.indexlist_complement:
            #    print rownames[myindex]
            within_cluster_num_elements = 0
            extra_cluster_num_elements = 0

            for i in range(len(thiscluster.indexlist)):

                this_column = variance_matrix[i,:]
                #print "this column in indexlist"
                #print this_column[thiscluster.indexlist]

                
                this_row = variance_matrix[:,i]
                #print "this row in indexlist"
                #print this_row[thiscluster.indexlist]
                thiscluster.intra_cluster_variance += sum(this_column[thiscluster.indexlist])
                thiscluster.extra_cluster_variance += sum(this_column[thiscluster.indexlist_complement])
                within_cluster_num_elements += len(thiscluster.indexlist)
                extra_cluster_num_elements += len(thiscluster.indexlist_complement)
                #within_cluster_num_elements += len(thiscluster.indexlist)

                #use all values not in this cluster
            #for j in range(len(thiscluster.indexlist_complement)):
            #    for k in range(len(thiscluster.indexlist_complement)):
            #        thiscluster.extra_cluster_variance += sum(variance_matrix[j,k])
            #        extra_cluster_num_elements +=1
                #extra_cluster_num_elements += len(thiscluster.indexlist_complement)
                
                
                #print "median var within cluster :"+str(median(this_row[thiscluster.indexlist]))
                #print "median var outside cluster :"+str(median(this_row[thiscluster.indexlist_complement]))
                

            
            thiscluster.intra_cluster_variance /= within_cluster_num_elements
            thiscluster.extra_cluster_variance /= extra_cluster_num_elements
            print "cluster "+str(clusternum)+", elements: "+str((within_cluster_num_elements))+",outside:"+str((extra_cluster_num_elements))+" within cluster average variance "+str(thiscluster.intra_cluster_variance)+" , outside of cluster average variance "+str(thiscluster.extra_cluster_variance)
            outfile.write("alter cluster"+str(clusternum)+", b="+str(thiscluster.intra_cluster_variance / thiscluster.extra_cluster_variance )+"\n")
            
            #outfile.write("color "+str(int(self.number-1))+", cluster"+str(self.number)+"\n"
    
    
    #finish session file
    outfile.write("cmd.show('cartoon'   ,'system')\n")
    if options.matrix != None:
        outfile.write("preset.b_factor_putty('system')"+"\n")  #,_self=cmd"+"\n")
        outfile.write("cmd.spectrum('b',selection=('all'),quiet=0)"+" \n")
        outfile.write("sele b_gt1, b > 1.0 \n")
        outfile.write("alter b_gt1, b = 1.0 \n")
        outfile.write("color bluewhite, b_gt1"+" \n")
    outfile.write("cmd.show('sticks','((byres (system))&(!(n;c,o,h|(n. n&!r. pro))))')" + " \n")
    if options.matrix != None:
        outfile.write("cmd.spectrum('b',selection=('all'),quiet=0)"+" \n")
        outfile.write("color bluewhite, b_gt1"+" \n")
    outfile.write("cmd.hide('(all and hydro)')\n")
    outfile.write("alter system, resi=str(int(resi)+"+str(options.begin)+")\n")
    outfile.write("\n")
    outfile.close()
    # Run community analysis
    #usage="%prog [-t traj1:traj2] [-x xvg_basedir] resfile [simulation numbers to use]  # where resfile is in the format <1-based-index> <aa type> <res num>"
    
    #communities_options=OptionParser(usage)
    #communities_options.add_option("-s", "--structure", default=None, type="string", help="pdb file")
    #communities_options.add_option("-f", "--filename", default=None, type="string", help="mutual information matrix")
    #communities_options.add_option("-o", "--outfile", default="communities", type="string", help="prefix for output files")
    #communities_options.add_option("-b", "--begin", default=0, type=int, help="first residue offset")
    #communities_options.add_option("-m", "--matrix", default=None, type="string", help="matrix for contacts")
    #communities_options.add_option("-e", "--edgelist", default="edgelist.dat",type="string", help="output edgelist data file")
    #communities_options.add_option("-c", "--cutoff", default=0.5, type="float", help="cutoff value for edges")
    #communities_options.add_option("-p", "--prefix",default="mutinf", type="string", help="prefix for output filenames")
    #communities_options.add_option("-i", "--iterations",default=None, type="int", help="number of iterations")
    #communities_options.add_option("-d", "--display_structure",default=None, type="string", help="full structure for pymol display")
    
    class communities_run_options:
        structure = None
        filename = None
        outfile = "communities"
        begin = 0
        matrix = None
        edgelist = "edgelist.dat"
        cutoff = 0.5
        prefix = "mutinf"
        iterations = None
        display_structure = None
        nocontacts = "no"
        user_communities = None

        def __init__(self):
            return
    
    
    ## Community analysis on mutinf matrix
    if(options.mutinf != "none"):
        communities_options = communities_run_options()
        communities_options.structure = options.hotatoms
        communities_options.filename = options.mutinf
        #communities_options.outfile = options.prefix + "_communities_mutinf"
        communities_options.begin = options.begin
        communities_options.matrix = options.contacts
        #communities_options.edgelist = options.prefix +"_communities_mutinf_edgelist.dat"
        communities_options.cutoff = options.cutoff
        #communities_options.prefix = options.prefix + "_mutinf"
        communities_options.iterations = options.iterations
        communities_options.display_structure = options.structure
        communities_options.nocontacts = options.nocontacts
        #communities_options.prefix = options.prefix
        communities_options.equal_weights = options.equal_weights
        communities_options.user_communities = options.user_communities
        # run community analysis from communities.py
        mylambda = lambda(x): -.5*log(1-exp(-2*x/3))
        run_communities(communities_options, mylambda)



    ## Community analysis on dist variance matrix
    if(options.matrix != "none"):
        communities2_options = communities_run_options()

        communities2_options.structure = options.hotatoms
        communities2_options.filename = options.matrix  #difference here
        communities2_options.outfile = options.prefix + "_communities_dist_variance"
        communities2_options.begin = options.begin
        communities2_options.matrix = options.contacts
        communities2_options.edgelist = options.prefix + "_communities_dist_variance_dist_variance_edgelist.dat"
        communities2_options.cutoff = 0.0
        communities2_options.prefix = options.prefix + "_dist_variance"
        communities2_options.iterations = options.iterations
        communities2_options.display_structure = options.structure
        communities2_options.nocontacts = options.nocontacts

        # run community analysis from communities.py
        mylambda = lambda(x): -.5*log(exp(-2*x/3))
        #run_communities(communities2_options, mylambda)
    
    
    

#END
        


