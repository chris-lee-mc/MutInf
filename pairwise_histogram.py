## Program to plot 2D histograms for torsions with highest mutual information given a pair of residues

import sys
sys.path.append('/home/cmcclend/mutinf/')
from input_output import *
from numpy import *
import matplotlib
import matplotlib.pyplot as plt

def xuniqueCombinations(items, n): # takes n distinct elements from the sequence, order is irrelevant
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
		    yield [items[i]]+cc


class Mutinf_Matrix_Chis_Bootstraps:
    mut_info_res_matrix = None
    rownames = []
    colnames = []
    old_bootstrap_sets = 0
    prefix = ""

    def __init__(self, prefix, nbins):
        self.old_bootstrap_sets = 0 #figure out # of bootstrap sets originally from output data
        self.prefix = prefix
        while (os.path.exists(self.prefix+"_bootstrap_"+str(self.old_bootstrap_sets)+"_mutinf_0diag.txt")):
           self.old_bootstrap_sets += 1
        if self.old_bootstrap_sets == 0:
           print "cannot find file"+self.prefix+"run1"+str(0)+"_mutinf_0diag.txt"
        print "old bootstrap sets: "+str(self.old_bootstrap_sets)
        ## READ MATRIX
        # self.prefix = "%s-nsims%d-structs%d-bin%d" % (os.path.basename(self.resfile_fn), self.num_sims, self.num_structs, int(self.binwidth))
        (test_matrix, self.rownames, self.colnames) = read_matrix_chis(self.prefix+"_bootstrap_0_mutinf_0diag.txt")
        self.rownames = self.colnames
        ##print test_matrix
        self.twoD_hist_boots = zeros((self.old_bootstrap_sets, test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3],nbins,nbins),float64)
        print "shape of twoD histograms over bootstrap sets: "+str(shape(self.twoD_hist_boots))
        for bootstrap in range(self.old_bootstrap_sets):
            twoD_temp = zeros((test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3], nbins, nbins),float64) #use as a temporary matrix to fill in the next line
            (self.twoD_hist_boots[bootstrap,:,:,:,:,:,:], self.rownames, self.colnames) = read_matrix_chis_2dhists(self.prefix+"_bootstrap_"+str(bootstrap)+"_2d_hists.txt", twoD_temp, nchi=6, nbins=nbins)
   
    def plot_pairwise_histogram(self,matrix,nbins,rownames,colnames,prefix):
       fig = plt.figure(figsize=(8,10),facecolor='w')
       ax1 = fig.add_axes([0.1,0.85,0.8,0.08],frameon=False)
       #dend = sch.dendrogram(links,orientation='top')
       ax1.set_xticks([])
       ax1.set_yticks([])
       mybins = arange(nbins)
       axmatrix = fig.add_axes([0.1,0.1,0.8,0.7])
       #idx = dend['leaves']
       #reordered = self.mymatrix[idx,:]
       #reordered = reordered[:,idx]
       im = axmatrix.matshow(matrix, aspect='auto', origin='upper', cmap=plt.cm.binary)
       #im.set_clim(0.0,max(matrix[:,:]))
       axmatrix.set_xticks(mybins)
       axmatrix.set_xticklabels(arange(0,360,360/nbins),rotation='vertical')
       axmatrix.set_yticks(mybins)
       axmatrix.set_yticklabels(arange(0,360,360/nbins),rotation='vertical')
       #for tick in axmatrix.xaxis.get_major_ticks():
           #tick.tick1On = False
           #tick.tick2On = False
       for label in axmatrix.xaxis.get_ticklabels():
           label.set_rotation(90)
           label.set_size(6)
       for label in axmatrix.yaxis.get_ticklabels():
           label.set_rotation(90)
           label.set_size(6)
       fig.suptitle('Pairwise Histogram',size=12.0)
       plt.plotname = str(prefix)+"_bootstrap_"+str(mybootstrap)+"_2d_hist_plot"
       plt.savefig(plt.plotname+".pdf",format="pdf")
       #os.system('ps2pdf'+' '+plt.plotname)
       #os.system('rm '+' '+plt.plotname+".ps")
       return fig

if __name__ == "__main__":
    usage="%prog resfile residue_i residue_j chi_i chi_j  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    #parser.add_option("-x", "--xvg_basedir", default=None, type="string", help="basedir to look for xvg files")
    #parser.add_option("-s", "--sigalpha", default=0.01, type="float", help="p-value threshold for statistical filtering, lower is stricter")
    parser.add_option("-s", "--num_structs", default=10001, type=int, help="number of snapshots per sim")
    parser.add_option("-w", "--binwidth", default=15.0, type="float", help="width of the bins in degrees")
    parser.add_option("-n", "--num_sims", default=None, type="int", help="number of simulations")
    #parser.add_option("-p", "--permutations", default=0, type="int", help="number of permutations for independent mutual information, for subtraction from total Mutual Information")
    #parser.add_option("-d", "--xvg_chidir", default = "/dihedrals/g_chi/", type ="string", help="subdirectory under xvg_basedir/run# where chi angles are stored")
    #parser.add_option("-a", "--adaptive", default = "yes", type ="string", help="adaptive partitioning (yes|no)")
    #parser.add_option("-b", "--backbone", default = "phipsichi", type = "string", help="chi: just sc  phipsi: just bb  phipsichi: bb + sc")
    parser.add_option("-o", "--bootstrap_set_size", default = None, type = "int", help="perform bootstrapping within this script; value is the size of the subsets to use")
    #parser.add_option("-i", "--skip", default = 1, type = "int", help="interval between snapshots to consider, in whatever units of time snapshots were output in") 
    #parser.add_option("-c", "--correct_formutinf_between_sims", default = "no", type="string", help="correct for excess mutual information between sims")
    #parser.add_option("-l", "--load_matrices_numstructs", default = 0, type = "int", help="if you want to load bootstrap matrices from a previous run, give # of structs per sim")
    #parser.add_option("--plot_2d_histograms", default = False, action = "store_true", help="makes 2d histograms for all pairs of dihedrals in the first bootstrap")
    #parser.add_option("-z", "--zoom_to_step", default = 0, type = "int", help="skips the first n snapshots in xvg files")
    #parser.add_option("-m","--max_num_chis", default = 99, type = "int", help="max number of sidechain chi angles per residue or ligand")
    #parser.add_option("-g","--gcc", default = 'gcc', type = "string", help="numpy distutils ccompiler to use. Recommended ones intelem or gcc")
    (options,args)=parser.parse_args()
    
    resfile_fn = args[0]
    residue_i = args[1]
    residue_j = args[2]
    chi_i = args[3]
    chi_j = args[4]
    prefix = "%s-nsims%d-structs%d-bin%d" % (resfile_fn, options.num_sims, options.num_structs, int(options.binwidth))
    print "prefix: " + str(prefix)
    nbins = int(360.0 / options.binwidth)
    #get number of bootstrap samples
    which_runs = []
    pair_runs_list = []
    for myruns in xuniqueCombinations(range(options.num_sims), options.bootstrap_set_size):
        which_runs.append(myruns)
    print which_runs
    #

    mymatrix = Mutinf_Matrix_Chis_Bootstraps(prefix, nbins=nbins)
    bootstrap_sets = mymatrix.old_bootstrap_sets
    Pij = mymatrix.twoD_hist_boots[:,residue_i, residue_j, chi_i, chi_j ]
    figlist = []

    for mybootstrap in range(bootstrap_sets):
	    figlist.append( mymatrix.plot_pairwise_histogram(Pij[mybootstrap],nbins,mymatrix.rownames,mymatrix.colnames,prefix))
            
	    
    #plt.show()
    
    

    
