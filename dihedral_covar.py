
from numpy import *
#import mdp
import re, os, sys, os.path, time, shelve
from optparse import OptionParser
from scipy import weave
from scipy import stats as stats
from scipy import special as special
from scipy import integrate as integrate
from scipy import misc as misc
from scipy.weave import converters
import time
import PDBlite, utils
from triplet import *
from constants import *
from input_output import *
from scipy.stats import gaussian_kde
from scipy.integrate import dblquad
from scipy.integrate import quad
from scipy.linalg.fblas import dger, dgemm

from dihedral_mutent import *

try:
       import MDAnalysis 
except: pass

set_printoptions(linewidth=120)

def cov2cor(cov):
    n=(shape(cov))[0]
    sigma=[0]*n
    cor=zeros((n,n),float64)
    for i in range(n):
        sigma[i]=sqrt(cov[i][i])
        for j in range(0, i+1):
            cor[i][j]=cor[j][i]=cov[i][j]/(sigma[i]*sigma[j] + SMALL*SMALL)
            pass
        pass
    return cor

#these two routines overwrite those in dihedral_mutent.py
def load_resfile(run_params, load_angles=True, all_angle_info=None):
    rp = run_params
    if rp.num_structs == None: rp.num_structs =  16777216 # # 1024 * 1024 * 16 #1500000
    sequential_num = 0
    resfile=open(rp.resfile_fn,'r')
    reslines=resfile.readlines()
    resfile.close()
    reslist = []
    for resline in reslines:
       if len(resline.strip()) == 0: continue
       xvg_resnum, res_name, res_numchain = resline.split()
       myexpr = re.compile(r"([0-9]+)([A-Z]*)")
       matches = myexpr.match(res_numchain)
       res_num = matches.group(1)
       if matches.group(2) != None:
              res_chain = matches.group(2)
       else:
              res_chain = ""
       if load_angles: 
              reslist.append(StressChis(res_name,res_num, res_chain, xvg_resnum, rp.xvg_basedir, rp.num_sims, rp.num_structs, rp.xvgorpdb, rp.binwidth, rp.sigalpha, rp.permutations, rp.phipsi, rp.backbone_only, rp.adaptive_partitioning, rp.which_runs, rp.pair_runs, bootstrap_choose = rp.bootstrap_choose, calc_variance=rp.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp.xvg_chidir, skip=rp.skip,skip_over_steps=rp.skip_over_steps,last_step=rp.last_step, calc_mutinf_between_sims=rp.calc_mutinf_between_sims,max_num_chis=rp.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp.pdbfile, xtcfile=rp.xtcfile, output_timeseries=rp.output_timeseries, lagtime_interval=rp.lagtime_interval, markov_samples=rp.markov_samples, num_convergence_points=rp.num_convergence_points ))
       else:  reslist.append(ResListEntry(res_name,res_num,res_chain))
       sequential_num += 1 
    return reslist

def load_data(run_params):
    load_resfile(run_params, load_angles=False) # make sure the resfile parses correctly (but don't load angle data yet)

    ### load trajectories from pdbs
    all_angle_info = None
    if run_params.xvgorpdb == "pdb":
       trajs = [PDBlite.PDBTrajectory(traj_fn) for traj_fn in run_params.traj_fns]
       traj_lens = array([traj.parse_len() for traj in trajs])
       run_params.num_structs = int(min(traj_lens)) # convert to python int, because otherwise it stays as a numpy int which weave has trouble interpreting
       run_params.num_res = trajs[0].parse_len()

       all_angle_info = AllAngleInfo(run_params)
       runs_to_load = set(array(run_params.which_runs).flatten())
       for sequential_sim_num, true_sim_num in zip(range(len(runs_to_load)), runs_to_load):
          all_angle_info.load_angles_from_traj(sequential_sim_num, trajs[true_sim_num-1], run_params, CACHE_TO_DISK)
       print "Shape of all angle matrix: ", all_angle_info.all_chis.shape

    print run_params
    print type(run_params.num_structs)

    ### load the residue list and angle info for those residues
    ### calculate intra-residue entropies and their variances
    reslist = load_resfile(run_params, load_angles=True, all_angle_info=all_angle_info)

    print "\n--Num residues: %d--" % (len(reslist))
    if (run_params.xvgorpdb): run_params.num_structs = int(max(reslist[0].numangles))

    return reslist



master_angles_matrix=None
test_reslist = None

class StressChis(ResidueChis):

   def __init__(self,myname,mynum,mychain,xvg_resnum,basedir,num_sims,max_angles,xvgorpdb,binwidth,sigalpha=1,
                permutations=0,phipsi=0,backbone_only=0,adaptive_partitioning=0,which_runs=None,pair_runs=None,bootstrap_choose=3,
                calc_variance=False, all_angle_info=None, xvg_chidir = "/dihedrals/g_chi/", skip=1, skip_over_steps=0, last_step=None, calc_mutinf_between_sims="yes", max_num_chis=99,
                sequential_res_num = 0, pdbfile = None, xtcfile = None, output_timeseries = "no", minmax=None, bailout_early = False, lagtime_interval = None, markov_samples = 250, num_convergence_points=1):
          
      global master_angles_matrix
      global test_reslist
      global xtc_coords 
      global last_good_numangles # last good value for number of dihedrals
      global NumChis, NumChis_Safe
      self.name = myname
      self.num = mynum
      self.chain = mychain
      self.xvg_basedir = basedir
      self.xvg_chidir = xvg_chidir
      self.xvg_resnum = xvg_resnum
      self.sequential_res_num = sequential_res_num
      self.backbone_only, self.phipsi = backbone_only, phipsi
      self.max_num_chis = max_num_chis
      self.markov_samples = markov_samples
      coarse_discretize = None
      split_main_side = None

      # we will look at mutual information convergence by taking linear subsets of the data instead of bootstraps, but use the bootstraps data structures and machinery. The averages over bootstraps then won't be meaningful
      # however the highest number bootstrap will contain the desired data -- this could be fixed later at the bottom of the code if desired
      # I also had to change some things in routines above that this code references in order to change numangles_bootstrap. We will essentially look at convergence by only looking at subsets of the data
      # in the weaves below, numangles will vary with 


      if(phipsi >= 0): 
             try:
                    self.nchi = self.get_num_chis(myname) * (1 - backbone_only) + phipsi * self.has_phipsi(myname)
             except:
                    NumChis = NumChis_Safe #don't use Ser/Thr hydroxyls for pdb trajectories
                    self.nchi = NumChis[myname] * (1 - backbone_only) + phipsi * self.has_phipsi(myname)
      elif(phipsi == -2):
             split_main_side = True
             if(self.chain == "S"):
                    self.nchi =  self.get_num_chis(myname)
             else:
                    self.nchi = 2 * self.has_phipsi(myname)
      elif(phipsi == -3):
             self.nchi = 3 #C-alpha x, y, z
      elif(phipsi == -4):
             print "doing analysis of stress data"
             self.nchi = 1 # just phi as a placeholder for a single variable
      else:             #coarse discretize phi/psi into 4 bins: alpha, beta, turn, other
             self.nchi = self.get_num_chis(myname) * (1 - backbone_only) + 1 * self.has_phipsi(myname)
             coarse_discretize = 1
             phipsi = 1
      if(xtcfile != None):
          self.nchi = 3 # x, y, z
      self.symmetry = ones((self.nchi),int16)
      self.numangles = zeros((num_sims),int32)
      self.num_sims = num_sims
      self.which_runs = array(which_runs)
      which_runs = self.which_runs
      #which_runs=array(self.which_runs)
      self.pair_runs = pair_runs
      self.permutations= permutations
      self.calc_mutinf_between_sims = calc_mutinf_between_sims
      if(bootstrap_choose == 0):
        bootstrap_choose = num_sims
      #print "bootstrap set size: "+str(bootstrap_choose)+"\n"
      #print "num_sims: "+str(num_sims)+"\n"
      #print self.which_runs
      #print "\n number of bootstrap sets: "+str(len(self.which_runs))+"\n"

      #check for free memory at least 15%
      #check_for_free_mem()
      
      #allocate stuff
      bootstrap_sets = self.which_runs.shape[0]

      #check num convergence points
      if num_convergence_points > 1:
             assert(num_convergence_points == bootstrap_sets)

      self.entropy =  zeros((bootstrap_sets,self.nchi), float64)
      self.entropy2 =  zeros((bootstrap_sets,self.nchi), float64) #entropy w/fewer bins
      self.entropy3 =  zeros((bootstrap_sets,self.nchi), float64) #entropy w/fewer bins
      self.entropy4 =  zeros((bootstrap_sets,self.nchi), float64) #entropy adaptive
      self.var_ent =  zeros((bootstrap_sets,self.nchi), float64)
      self.numangles_bootstrap = zeros((bootstrap_sets),int32)
      print "\n#### Residue: "+self.name+" "+self.num+" "+self.chain+" torsions: "+str(self.nchi), utils.flush()
      binwidth = float(binwidth)
      bins = arange(0,360, binwidth) #  bin edges global variable
      nbins=len(bins) # number of bins
      nbins_cor = int(nbins * FEWER_COR_BTW_BINS);
      self.nbins=nbins
      self.nbins_cor=nbins_cor
      sqrt_num_sims=sqrt(num_sims)
      self.chi_pop_hist=zeros((bootstrap_sets, self.nchi,nbins),float64)
      self.chi_counts=zeros((bootstrap_sets, self.nchi, nbins), float64) #since these can be weighted in advanced sampling
      #self.chi_var_pop=zeros((bootstrap_sets, self.nchi,nbins),float64)
      self.chi_pop_hist_sequential=zeros((num_sims, self.nchi, nbins_cor), float64)
      num_histogram_sizes_to_try = 2  # we could try more and pick the optimal size
      self.chi_counts_sequential=zeros((num_sims, self.nchi, nbins_cor), float64) #half bin size
      self.chi_counts_sequential_varying_bin_size=zeros((num_histogram_sizes_to_try, num_sims, self.nchi, int(nbins*(num_histogram_sizes_to_try/2)) ), float64) #varying bin size
      self.angles_input = zeros((self.nchi,num_sims,max_angles),float64)         # the dihedral angles, with a bigger array than will be needed later

      #self.sorted_angles = zeros((self.nchi,num_sims,max_angles),float64) # the dihedral angles sorted
      self.ent_hist_left_breaks = zeros((self.nchi, nbins * MULT_1D_BINS + 1),float64)
      self.adaptive_hist_left_breaks = zeros((bootstrap_sets, nbins + 1),float64) #nbins plus one to define the right side of the last bin
      self.adaptive_hist_left_breaks_sequential = zeros(( num_sims, nbins_cor + 1 ),float64)  #nbins_cor plus one to define the right side of the last bin
      self.adaptive_hist_binwidths = zeros((bootstrap_sets, nbins ),float64) 
      self.adaptive_hist_binwidths_sequential = zeros(( num_sims, nbins_cor ),float64)
      self.ent_hist_binwidths = zeros((bootstrap_sets, self.nchi, nbins * MULT_1D_BINS),float64)
      self.ent_from_sum_log_nn_dists = zeros((bootstrap_sets, self.nchi, MAX_NEAREST_NEIGHBORS),float64)
      self.minmax = zeros((2,self.nchi))
      self.minmax[1,:] += 1 #to avoid zero in divide in expand_contract angles



      if(phipsi >= 0): 
             self.nchi = self.get_num_chis(myname) * (1 - backbone_only) + phipsi * self.has_phipsi(myname)
      elif(phipsi == -2):
             split_main_side = True
             if(self.chain == "S"):
                    self.nchi =  self.get_num_chis(myname)
             else:
                    self.nchi = 2 * self.has_phipsi(myname)
      elif(phipsi == -3):
             self.nchi = 3 #C-alpha x, y, z
      elif(phipsi == -4):
             print "doing analysis of stress data"
             self.nchi = 1 # just phi as a placeholder for a single variable
      else:             #coarse discretize phi/psi into 4 bins: alpha, beta, turn, other
             self.nchi = self.get_num_chis(myname) * (1 - backbone_only) + 1 * self.has_phipsi(myname)
             coarse_discretize = 1
             phipsi = 1
      if(xtcfile != None):
          self.nchi = 3 # x, y, z
      self.symmetry = ones((self.nchi),int16)
      self.numangles = zeros((num_sims),int32)
      self.num_sims = num_sims
      self.which_runs = array(which_runs)
      which_runs = self.which_runs
      #which_runs=array(self.which_runs)
      self.pair_runs = pair_runs
      self.permutations= permutations
      self.calc_mutinf_between_sims = calc_mutinf_between_sims
      if(bootstrap_choose == 0):
        bootstrap_choose = num_sims
      #print "bootstrap set size: "+str(bootstrap_choose)+"\n"
      #print "num_sims: "+str(num_sims)+"\n"
      #print self.which_runs
      #print "\n number of bootstrap sets: "+str(len(self.which_runs))+"\n"

      #check for free memory at least 15%
      #check_for_free_mem()
      
      #allocate stuff
      bootstrap_sets = self.which_runs.shape[0]

      #check num convergence points
      if num_convergence_points > 1:
             assert(num_convergence_points == bootstrap_sets)

      
      self.numangles_bootstrap = zeros((bootstrap_sets),int32)
      print "\n#### Residue: "+self.name+" "+self.num+" "+self.chain+" torsions: "+str(self.nchi), utils.flush()

      if(xvgorpdb == "xvg"):
         self._load_xvg_data(basedir, num_sims, max_angles, xvg_chidir, skip,skip_over_steps,last_step, coarse_discretize, split_main_side)
      if(xvgorpdb == "pdb"):
         self._load_pdb_data(all_angle_info, max_angles)
      if(xvgorpdb == "xtc"):
         self.load_xtc_data(basedir, num_sims, max_angles, xvg_chidir, skip, skip_over_steps, pdbfile, xtcfile)

      #resize angles array to get rid of trailing zeros, use minimum number
      print "weights"
      print self.weights

      #print "resizing angles array, and creating arrays for adaptive partitioning" 
      min_num_angles = int(min(self.numangles))
      max_angles = int(min_num_angles)
      if(min_num_angles > 0):
             last_good_numangles = min_num_angles
      self.angles = zeros((self.nchi, num_sims, min_num_angles))
      angles_autocorrelation = zeros((self.nchi, bootstrap_sets, min_num_angles), float64)
      #bins_autocorrelation =   zeros((self.nchi, bootstrap_sets, min_num_angles), float64)
      self.boot_sorted_angles = zeros((self.nchi,bootstrap_sets,bootstrap_choose*max_angles),float64)
      self.boot_ranked_angles = zeros((self.nchi,bootstrap_sets,bootstrap_choose*max_angles),int32) 
      self.boot_weights = zeros((bootstrap_sets,bootstrap_choose*max_angles),float64) 
      self.rank_order_angles = zeros((self.nchi,num_sims,max_angles),int32) # the dihedral angles
                                                                        # rank-ordered with respect to all sims together
      self.rank_order_angles_sequential = zeros((self.nchi,num_sims,max_angles),int32) # the dihedral angles
                                                                        # rank-ordered for each sim separately
                                                                        # for mutinf between sims
      #max_num_angles = int(max(self.numangles))
      max_num_angles = int(min(self.numangles)) #to avoid bugs
      

      #counts_marginal=zeros((bootstrap_sets,self.nchi,nbins),float32) # normalized number of counts per bin, 


      #print "initialized angles_new array"
      self.numangles[:] = min(self.numangles)
      #print "new numangles"
      #print self.numangles
      for mychi in range(self.nchi):
             for num_sim in range(num_sims):
                    self.angles[mychi,num_sim,:min_num_angles] = self.angles_input[mychi,num_sim,:min_num_angles]
      #print "done copying angles over"
      del self.angles_input #clear up memory space
      #print self.angles
      
      
      

      if ((self.name in ("GLY", "ALA")) and (phipsi == 0 or phipsi == -2)): 
           # First prepare chi pop hist and chi pop hist sequential needed for mutual information -- just dump everything into one bin, 
           # giving entropy zero, which should also give MI zero, but serve as a placeholder in the mutual information matrix
           if(last_good_numangles > 0):  
                  self.numangles[:] = last_good_numangles # dummy value  
           self.numangles_bootstrap[:] = bootstrap_choose *  int(min(self.numangles))
           
                  
           
           return   #if the side chains don't have torsion angles, drop out

      self.sorted_angles = zeros((self.nchi, sum(self.numangles)),float64)
      
      if(xvgorpdb == "xvg" or (xvgorpdb == "pdb" and phipsi != -3)): #if not using C-alphas from pdb
             self.correct_and_shift_angles(num_sims,bootstrap_sets,bootstrap_choose, coarse_discretize)
      elif(xvgorpdb == "xtc" or (xvgorpdb == "pdb" and phipsi == -3)) : #if using xtc cartesians or pdb C-alphas 
             self.correct_and_shift_carts(num_sims,bootstrap_sets,bootstrap_choose, num_convergence_points)
             
      if(minmax == None):
             print "getting min/max values"
             mymin = zeros(self.nchi)
             mymax = zeros(self.nchi)
             for mychi in range(self.nchi):
                    mymin[mychi] =  min((self.angles[mychi,:,:min(self.numangles)]).flatten())
                    mymax[mychi] =  max((self.angles[mychi,:,:min(self.numangles)]).flatten())
             print "__init__ mymin: "
             print mymin
             print "__init__ mymax: "
             print mymax
             self.minmax[0, :] = mymin
             self.minmax[1, :] = mymax
             for mychi in range(self.nchi):
                    if(self.minmax[1,mychi] - self.minmax[0,mychi] <= 0):
                           self.minmax[1,mychi] = self.minmax[0,mychi] + 1
             print self.minmax
      else:
             self.minmax = minmax
             self.expand_contract_data(num_sims,bootstrap_sets,bootstrap_choose)
      
      if(master_angles_matrix == None):
          master_angles_matrix=zeros((run_params.num_sims, len(test_reslist), min_num_angles),float64) #nchi is 1 for stress analysis

      
      master_angles_matrix[:,sequential_res_num,:] = self.angles[0,:,:]
      print "master angles matrix: "
      print master_angles_matrix[0]
      del self.angles #cleanup
      del self.rank_order_angles
      del self.rank_order_angles_sequential
      del self.sorted_angles
      del self.boot_sorted_angles
      del self.boot_ranked_angles
      








##########################################################################################################    
#===================================================
#READ INPUT ARGUMENTS
#===================================================
##########################################################################################################


if __name__ == "__main__":
    try:
       import run_profile
       run_profile.fix_args()
    except: pass
    
    usage="%prog [-t traj1:traj2] [-x xvg_basedir] resfile [simulation numbers to use]  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    parser.add_option("-t", "--traj_fns", default=None, type="string", help="filenames to load PDB trajectories from; colon separated (e.g. fn1:fn2)")
    parser.add_option("-x", "--xvg_basedir", default=None, type="string", help="basedir to look for xvg files")
    parser.add_option("-s", "--sigalpha", default=0.01, type="float", help="p-value threshold for statistical filtering, lower is stricter")
    parser.add_option("-w", "--binwidth", default=15.0, type="float", help="width of the bins in degrees")
    parser.add_option("-n", "--num_sims", default=None, type="int", help="number of simulations")
    parser.add_option("-p", "--permutations", default=0, type="int", help="number of permutations for independent mutual information, for subtraction from total Mutual Information")
    parser.add_option("-d", "--xvg_chidir", default = "/dihedrals/g_chi/", type ="string", help="subdirectory under xvg_basedir/run# where chi angles are stored")
    parser.add_option("-a", "--adaptive", default = "yes", type ="string", help="adaptive partitioning (yes|no)")
    parser.add_option("-b", "--backbone", default = "phipsichi", type = "string", help="chi: just sc  phipsi: just bb  phipsichi: bb + sc")
    parser.add_option("-o", "--bootstrap_set_size", default = None, type = "int", help="perform bootstrapping within this script; value is the size of the subsets to use")
    parser.add_option("-i", "--skip", default = 1, type = "int", help="interval between snapshots to consider, in whatever units of time snapshots were output in") 
    parser.add_option("-c", "--correct_formutinf_between_sims", default = "no", type="string", help="correct for excess mutual information between sims")
    parser.add_option("--load_matrices_numstructs", default = 0, type = "int", help="if you want to load bootstrap matrices from a previous run, give # of structs per sim (yes|no)")
    parser.add_option("-l", "--last_step", default=None, type= "int", help="last step to read from input files, useful for convergence analysis")
    parser.add_option("--plot_2d_histograms", default = False, action = "store_true", help="makes 2d histograms for all pairs of dihedrals in the first bootstrap")
    parser.add_option("-z", "--zoom_to_step", default = 0, type = "int", help="skips the first n snapshots in xvg files")
    parser.add_option("-M","--markov_samples", default = 0, type = "int", help="markov state model samples to use for independent distribution")
    parser.add_option("-N","--max_num_lagtimes", default = 5000, type =	"int", help="maximum number of lagtimes for markov model")
    parser.add_option("-m","--max_num_chis", default = 99, type = "int", help="max number of sidechain chi angles per residue or ligand")
    parser.add_option("-f","--pdbfile", default = None, type = "string", help="pdb structure file for additional 3-coord cartesian per residue")
    parser.add_option("-q","--xtcfile", default = None, type = "string", help="gromacs xtc prefix in 'run' subdirectories for additional 3-coord cartesian per residue")
    parser.add_option("-g","--gcc", default = 'intelem', type = "string", help="numpy distutils ccompiler to use. Recommended ones intelem or gcc")
    parser.add_option("-e","--output_timeseries", default = "yes", type = "string", help="output corrected dihedral timeseries (requires more memory) yes|no ")  #default yes for covariance analysis
    parser.add_option("-y","--symmetry_number", default = 1, type = int, help="number of identical subunits in homo-oligomer for symmetrizing matrix")
    parser.add_option("-L","--lagtime_interval", default = None, type=int, help="base snapshot interval to use for lagtimes in Markov model of bin transitions")
    parser.add_option("-j","--offset", default = 0,type=int, help="offset for mutinf of (i,j) at (t, t - offset)")
    parser.add_option("--output_independent",default = 0, type=int, help="set equal to 1 to output independent mutinf values from markov model or multinomial distribution")
    parser.add_option("-C","--num_convergence_points", default = 0, type=int, help="for -n == -o , use this many subsets of the data to look at convergence statistics")   
    parser.add_option("-T","--triplet", default=None, type="string", help="wheter to perform triplet mutual information or not")
    (options,args)=parser.parse_args()
    mycompiler = options.gcc
    if len(filter(lambda x: x==None, (options.traj_fns, options.xvg_basedir))) != 1:
        parser.error("ERROR exactly one of --traj_fns or --xvg_basedir must be specified")

    print "COMMANDS: ", " ".join(sys.argv)
    
    
    # Initialize
    offset = options.offset
    resfile_fn = args[0]
    adaptive_partitioning = (options.adaptive == "yes")
    phipsi = 0
    backbone_only = 0
    if options.backbone == "calpha" or options.backbone == "calphas":
           if options.traj_fns != None:
                  phipsi = -3
                  backbone_only = 1
           else:
                  options.backbone = "phipsi"
    
    if options.backbone == "phipsichi":
        phipsi = 2;
        backbone_only =0
    if options.backbone == "phipsi":
        phipsi = 2;
        backbone_only = 1
    if options.backbone == "stress":
        phipsi = -4;
        NumChis["GLY"] = 1
        NumChis["ALA"] = 1
        backbone_only = 1
        print "performing stress analysis"
    if options.backbone == "coarse_phipsi":
        phipsi = -1
        backbone_only = 1
        print "overriding binwidth, using four bins for coarse discretized backbone"
        options.binwidth = 90.0
    if options.backbone == "split_main_side":
        print "treating backbone and sidechain separately according to residue list"
        phipsi = -2
        backbone_only = 0

    #phipsi = options.backbone.find("phipsi")!=-1
    #backbone_only = options.backbone.find("chi")==-1
    bins=arange(-180,180,options.binwidth) #Compute bin edges
    nbins = len(bins)

    if options.traj_fns != None:
       xvgorpdb = "pdb"
       traj_fns = options.traj_fns.split(":")
       num_sims = len(traj_fns)
    else:
        if(options.xtcfile != None):
            xvgorpdb = "xtc"
            num_sims = options.num_sims
            traj_fns = None
        else:
            xvgorpdb = "xvg"
            traj_fns = None
            num_sims = options.num_sims
 
    print "num_sims:"+str(num_sims)+"\n"

    #if(len(args) > 1):
    #    which_runs = map(int, args[1:])
    #    num_sims = len(which_runs)
    #else:
    assert(num_sims != None)
    #which_runs = range(num_sims)
    #which_runs = None

    if options.bootstrap_set_size == None:
        options.bootstrap_set_size = num_sims

    which_runs = []
    pair_runs_list = []
    for myruns in xuniqueCombinations(range(num_sims), options.bootstrap_set_size):
        which_runs.append(myruns)
    print which_runs

    if options.num_convergence_points > 1: #create bootstrap samples for number of convergence points, this will also set variables like bootstrap_sets  
           print "looking at mutual information convergence using "+str(options.num_convergence_points)+" convergence points"
           options.bootstrap_set_size = 1  #override options
           for convergence_point in range(options.num_convergence_points - 1):
                  which_runs.append(which_runs[0])


    NUM_LAGTIMES=options.max_num_lagtimes #overwrite global var
    OUTPUT_INDEPENDENT_MUTINF_VALUES = options.output_independent

    bootstrap_pair_runs_list = []
    for bootstrap in range((array(which_runs)).shape[0]):
        pair_runs_list = []
        for myruns2 in xcombinations(which_runs[bootstrap], 2):
            pair_runs_list.append(myruns2)
        bootstrap_pair_runs_list.append(pair_runs_list)

    pair_runs_array = array(bootstrap_pair_runs_list,int16)
    print pair_runs_array 

    
    

    #set num_structs = options.load_matrices in case we don't actually load any data, just want to get the filenames right for the bootstrap matrices

    if (options.load_matrices_numstructs > 0): num_structs = options.load_matrices_numstructs
    else: num_structs = None

    run_params = RunParameters(resfile_fn=resfile_fn, adaptive_partitioning=adaptive_partitioning, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins,
      bootstrap_set_size=options.bootstrap_set_size, sigalpha=options.sigalpha, permutations=options.permutations, num_sims=num_sims, num_structs=num_structs,
      binwidth=options.binwidth, bins=bins, which_runs=which_runs, xvgorpdb=xvgorpdb, traj_fns=traj_fns, xvg_basedir=options.xvg_basedir, calc_variance=False, xvg_chidir=options.xvg_chidir,bootstrap_choose=options.bootstrap_set_size,pair_runs=pair_runs_array,skip=options.skip,skip_over_steps=options.zoom_to_step,last_step=options.last_step,calc_mutinf_between_sims=options.correct_formutinf_between_sims,load_matrices_numstructs=options.load_matrices_numstructs,plot_2d_histograms=options.plot_2d_histograms,max_num_chis=options.max_num_chis,pdbfile=options.pdbfile, xtcfile=options.xtcfile, output_timeseries=options.output_timeseries,lagtime_interval=options.lagtime_interval, markov_samples=options.markov_samples, num_convergence_points=options.num_convergence_points)


    print run_params

    
    
    #====================================================
    #DO ANALYSIS
    #===================================================

    print "Calculating Entropy and Mutual Information"
    #must initialize global test_reslist so that we know how many residues there are 
    test_reslist = load_resfile(run_params, load_angles=False, all_angle_info=None) #see how many residues there are
    

    independent_mutinf_thisdof = zeros((run_params.permutations,1),float32)
    timer = utils.Timer()


    ### load angle data, calculate entropies and mutual informations between residues (and error-propagated variances)
    if run_params.bootstrap_set_size == None:
        if run_params.load_matrices_numstructs == 0: 
               reslist = load_data(run_params)
               print "TIME to load trajectories & calculate intra-residue entropies: ", timer
               timer=utils.Timer()

               prefix = run_params.get_logfile_prefix()
               ##========================================================
               ## Output Timeseries Data in a big matrix
               ##========================================================
               name_num_list = make_name_num_list(reslist)
               if(xvgorpdb == "xtc"):
                      print "calculating distance matrices and variance:"
                      output_distance_matrix_variances(len(run_params.which_runs),run_params.bootstrap_set_size,run_params.which_runs, reslist[0].numangles, reslist[0].numangles_bootstrap, name_num_list) #uses global xtc_coords for data
               if run_params.output_timeseries == "yes":
                      timeseries_chis_matrix = output_timeseries_chis(prefix+"_timeseries",reslist,name_num_list,run_params.num_sims)
                      print "TIME to output timeseries data: ", timer
                      timer=utils.Timer()
               
 
               print "TIME to calculate pair stats: ", timer
               timer=utils.Timer()

        prefix = run_params.get_logfile_prefix() + "_sims" + ",".join(map(str, sorted(which_runs)))
    else:
        runs_superset, set_size = run_params.which_runs, run_params.bootstrap_set_size
        if set_size > len(runs_superset[0]) or len(runs_superset) < 1:
            print "FATAL ERROR: invalid values for bootstrap set size '%d' from runs '%s'" % (set_size, runs_superset)
            sys.exit(1)

        run_params.calc_variance = False
        print "\n----- STARTING BOOTSTRAP RUNS: %s -----" % run_params
        if run_params.load_matrices_numstructs == 0:
            reslist = load_data(run_params)
        print "TIME to load trajectories & calculate intra-residue entropies: ", timer
        timer=utils.Timer()
        
        

        ##========================================================
        ## Output Timeseries Data in a big matrix
        ##========================================================
        prefix = run_params.get_logfile_prefix()

        if run_params.load_matrices_numstructs == 0: 
           
           name_num_list = make_name_num_list(reslist)

           if(xvgorpdb == "xtc" and run_params.load_matrices_numstructs == 0 ):
                         name_num_list = make_name_num_list(reslist)
                         print "calculating distance matrices and variance:"
                         output_distance_matrix_variances(len(run_params.which_runs),run_params.bootstrap_set_size,run_params.which_runs, reslist[0].numangles, reslist[0].numangles_bootstrap, name_num_list) #uses global xtc_coords for data, len(which_runs) gives number of bootstrap_sets
           #if run_params.output_timeseries == "yes":
           #       timeseries_chis_matrix = output_timeseries_chis(prefix+"_timeseries",reslist,name_num_list,run_params.num_sims)
           #       print "TIME to output timeseries data: ", timer
           #       timer=utils.Timer()

           #if run_params.load_matrices_numstructs == 0:
               #mut_info_res_matrix, uncorrected_mutinf_thisdof, corrections_mutinf_thisdof, mut_info_uncert_matrix, mut_info_res_matrix_different_sims, dKLtot_dresi_dresj_matrix, twoD_hist_boot_avg, twoD_hist_boots, twoD_hist_ind_boots, mut_info_norm_res_matrix = calc_pair_stats(reslist, run_params)
              
           numangles_sum=sum(reslist[0].numangles)
           min_numangles=min(reslist[0].numangles)    
           
           #### note self.angles = zeros((self.nchi,num_sims,min_numangles),float64)   
           
           master_cov = zeros((run_params.num_sims, len(reslist), len(reslist)), float64)
           
           #populate master matrix
           #for myindex, res in zip(range(len(reslist)), reslist):
           #    master_angles_matrix[:,myindex,:] = res.angles[0,:,:]
           
           def mean_cov(X):
               n,p = X.shape
               m = X.mean(axis=0)
               # covariance matrix with correction for rounding error
               # S = (cx'*cx - (scx'*scx/n))/(n-1)
               # Am Stat 1983, vol 37: 242-247.
               cx = X - m
               cxT_a = zeros((p-1,n),float64)
               cxT_b = zeros((p-1,n),float64)
               cxT_a[:,:] = (cx.T)[0:p-1,:]
               cxT_b[:,:] = (cx.T)[1:p,:]
               scx = cx.sum(axis=0)
               scx_op = dger(-1.0/n,scx,scx)
               S = dgemm(1.0, cx.T, cx.T, beta=1.0,
                       c=scx_op, trans_a=0, trans_b=1, overwrite_c=1)
               #S = dgemm(1.0, cxT_a, cxT_b, beta=1.0,
               #        c=scx_op, trans_a=0, trans_b=1, overwrite_c=1)
               S[:] *= 1.0/(n-1)
               return m,S.T
           
           master_cov_matrix = zeros((num_sims, len(reslist),len(reslist)),float64)
           
           print "shape of master_angles_matrix:"
           print shape(master_angles_matrix)
           print "shape of transpose of master angles matrix:"
           print shape(transpose(master_angles_matrix))
           
           
           for mysim in range(num_sims):
               mymean, master_cov_matrix[mysim,:,:] = mean_cov(transpose(master_angles_matrix[mysim]))
                                                               
           #print mut_info_res_matrix
           print "TIME to calculate covar matrix: ", timer
           timer=utils.Timer()
                                                               
                                                               
                                                               
           # create a master matrix
           #matrix_list += [calc_pair_stats(reslist, run_params)[0]] # add the entropy/mut_inf matrix to the list of matrices
           #bootstraps_mut_inf_res_matrix = zeros(list(matrix_list[0].shape) + [len(matrix_list)], float32)
           #for i in range(len(matrix_list)): bootstraps_mut_inf_res_matrix[:,:,i] = matrix_list[i]

           #mut_info_res_matrix, mut_info_uncert_matrix = bootstraps_mut_inf_res_matrix.mean(axis=2), bootstraps_mut_inf_res_matrix.std(axis=2)
           #prefix = run_params.get_logfile_prefix() + "_sims%s_choose%d" % (",".join(map(str, sorted(runs_superset))), set_size)
        
    ### output results to disk

    mycov_avg = average(master_cov_matrix,axis=0)
    mycor_avg = cov2cor(mycov_avg)
    
    output_diag(prefix+"_avg_msf.txt",mycov_avg,name_num_list)
    output_matrix(prefix+"_avg_cov.txt",mycov_avg, name_num_list,name_num_list,zero_diag=False) 
    output_matrix(prefix+"_avg_cov_0diag.txt",mycov_avg, name_num_list,name_num_list,zero_diag=True) 
    output_matrix(prefix+"_avg_corr.txt",mycor_avg, name_num_list,name_num_list,zero_diag=True) 
    output_matrix(prefix+"_avg_corr_0diag.txt",mycor_avg, name_num_list,name_num_list,zero_diag=True) 

    print mycor_avg
    for i in range(mycor_avg.shape[0]):
                  mycor_avg[i][i] = 0
    
    mycor_avg_flat = mycor_avg.flatten()
    print mycor_avg_flat
    mystd = std(mycor_avg_flat[mycor_avg_flat != 0.0])
    print "standard dev: "+str(mystd)
    #or i in range(mycor_avg_flat.shape[0]):
    #      if mycor_avg_flat[i] > 0:
    #       if mycor_avg_flat/mystd < 1.644853626951 :
    #                mycor_avg_flat[i] = 0 
                     
    mycor_avg_flat[abs(mycor_avg_flat + SMALL*SMALL)/(mystd + SMALL*SMALL) < 1.644853626951   ] = 0 # 1.645 sigma -> 90% of values zeroed, 3sigma -> ~99.8 % of values zeroed
    
    mycor_avg_sigma = reshape(mycor_avg_flat, shape(mycor_avg))
    output_matrix(prefix+"_avg_corr_1_6sigma.txt",mycor_avg_sigma, name_num_list,name_num_list,zero_diag=True) 

    mycor_avg_flat[abs(mycor_avg_flat + SMALL*SMALL)/(mystd + SMALL*SMALL) < 2.0   ] = 0 # 1.645 sigma -> 90% of values zeroed, 3sigma -> ~99.8 % of values zeroed
    
    mycor_avg_sigma = reshape(mycor_avg_flat, shape(mycor_avg))
    output_matrix(prefix+"_avg_corr_2sigma.txt",mycor_avg_sigma, name_num_list,name_num_list,zero_diag=True) 


    mycor_avg_flat[abs(mycor_avg_flat + SMALL*SMALL)/(mystd + SMALL*SMALL) < 2.575829303549  ] = 0 # 2.5 sigma -> 99% of values zeroed, 3sigma -> ~99.8 % of values zeroed
    
    mycor_avg_2sigma = reshape(mycor_avg_flat, shape(mycor_avg))
    output_matrix(prefix+"_avg_corr_2_6sigma.txt",mycor_avg_2sigma, name_num_list,name_num_list,zero_diag=True) 

    #for mysim in range(num_sims):
    #       mycor = cov2cor(master_cov_matrix[mysim,:,:])
    #       
    #       output_matrix(prefix+"_sim_"+str(mysim)+"_covar.txt",master_cov_matrix[mysim,:,:],name_num_list,name_num_list,zero_diag=True)    
    #       output_matrix(prefix+"_sim_"+str(mysim)+"_corr.txt",cov2cor(master_cov_matrix[mysim,:,:]),name_num_list,name_num_list,zero_diag=True)    

    #END

