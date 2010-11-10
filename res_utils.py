#!/usr/bin/python
from numpy import zeros, ones, arange, outer, inner, float32, float64, int16, sqrt, floor, mean, std, log, add, average, nonzero, cov, transpose, resize, array, reshape, repeat, transpose, max, array2string, sum, histogram2d, logical_and, random, sort, searchsorted, swapaxes, shape, int8, int32
import re, os, sys, os.path, time, shelve, commands
from optparse import OptionParser
from scipy import weave
from scipy import stats as stats
from scipy import special as special
from scipy import integrate as integrate
from scipy.weave import converters
import PDBlite

SMALL = 0.000000001
PI = 3.141592653589793238462643383279502884197
TWOPI = 2 * PI
ROOTPI = 1.772453850905516027298167483341
GAMMA_EULER = 0.57721566490153286060651209
K_NEAREST_NEIGHBOR = 1
CACHE_TO_DISK = True
CORRECT_FOR_SYMMETRY = 1
VERBOSE = 1
OFF_DIAG = 1 #set to "1" to look at off-diagonal terms
OUTPUT_2D_PROB_PLOTS = False # output 2d probability plots for all pairs of dihedral angles (Pij, PiPj, and Pij/PiPj); slow and writes *lots* of files
OUTPUT_DIAG = 1  #set to "1" to output flat files for diagonal matrix elements
EACH_BOOTSTRAP_MATRIX = 1 #set to "1" to output flat files for each bootstrapped matrix
MULT_1D_BINS = 3 #set to "3" or more for MULT_1D_BINS times the number of 1-D bins for 1-D entropy only; single and double already done by defualt
SAFE_1D_BINS = 7 # approximately 2PI bins so that continuous h(x) = H(x, discrete) + log (binwidth)
FEWER_COR_BTW_BINS = 0.5 #set to "0.5" for half as many 2-D bins for 2-D correlation between tors in different sims, nbins must be an even number for this to be appropriate!
MAX_NEAREST_NEIGHBORS = 1 # up to k nearest neighbors
SPEED = "vfast" # set to "vfast" to use the vfast_cov(); "fast" to use fast_cov(); "*" to use cov()



# Stores the parameters for an Entropy/Mutinf calculation run
# e.g. Access with runparams.permutations
class RunParameters(object):
    def __init__(self, **kwds):
        self.__dict__ = kwds
        # print self.__dict__

    # allow get access e.g. o.myattribute
    def __getattr__(self, name):
        try: return self.__dict__[name]
        except KeyError:
            print "FATAL ERROR: can't find key '%s' in RunParamaters object" % name
            sys.exit(1)

    # allow set access e.g. o.myattribute = 2
    def __setattr_(self, name, val):
        self.__dict__[name] = val

    def __str__(self):
        d = self.__dict__
        return "Num sims = %d; Runs being used = %s; Binwidth = %.1f; Permutations = %d; Num Structs = %s" % \
               (d['num_sims'], " ".join(map(str, d['which_runs'])), d['binwidth'], d['permutations'], d['num_structs'])

    def get_logfile_prefix(self):
        return "%s-nsims%d-structs%d-bin%d" % (os.path.basename(self.resfile_fn), self.num_sims, self.num_structs, int(self.binwidth))


class ResListEntry:
    name = 'XXX'
    num = 0
    chain = ' '
    def __init__(self,myname,mynum):
        self.name = myname
        self.num = mynum

def load_resfile(run_params, load_angles=True, all_angle_info=None):
    rp = run_params
    if rp.num_structs == None: rp.num_structs = 15000

    resfile=open(rp.resfile_fn,'r')
    reslines=resfile.readlines()
    resfile.close()
    reslist = []
    sequential_num = 0 # not used, right now...
    for resline in reslines:
       if len(resline.strip()) == 0: continue
       dummy, res_name, res_num = resline.split()
       if load_angles: reslist.append(ResidueChis(res_name,res_num, sequential_num, rp.xvg_basedir, rp.num_sims, rp.num_structs, rp.xvgorpdb, rp.binwidth, rp.sigma,
                                                  rp.permutations, rp.phipsi, rp.backbone_only, rp.adaptive_partitioning, rp.which_runs, rp.pair_runs, bootstrap_choose = rp.bootstrap_choose, calc_variance=rp.calc_variance,
                                                  all_angle_info=all_angle_info, xvg_chidir=rp.xvg_chidir, skip=rp.skip, calc_mutinf_between_sims=rp.calc_mutinf_between_sims))
       else: reslist.append(ResListEntry(res_name,int(res_num)))
           
       sequential_num += 1
       
    return reslist

# Load angle data and calculate intra-residue entropies
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
          all_angle_info.load_angles_from_traj(sequential_sim_num, trajs[true_sim_num-1], CACHE_TO_DISK)
       print "Shape of all angle matrix: ", all_angle_info.all_chis.shape

    print run_params
    print type(run_params.num_structs)

    ### load the residue list and angle info for those residues
    ### calculate intra-residue entropies and their variances
    reslist = load_resfile(run_params, load_angles=True, all_angle_info=all_angle_info)

    print "\n--Num residues: %d--" % (len(reslist))
    if (run_params.xvgorpdb): run_params.num_structs = int(max(reslist[0].numangles))

    return reslist


class AllAngleInfo:
   CACHE_FN = os.path.expanduser("~/dihed_mutent.cache")

   def __init__(self, run_params):
      rp = run_params
      self.all_chis = None
      self.num_sims, self.num_structs, self.num_res = rp.num_sims, rp.num_structs, rp.num_res
      self.calc_variance = False
   # load dihedrals from a pdb trajectory
   def load_angles_from_traj(self, sequential_sim_num, pdb_traj, cache_to_disk=True):
      sim_chis, found_in_cache = None, False

      # check if the simulations is cached
      if cache_to_disk:
          key = pdb_traj.traj_fn
          if os.path.exists(self.CACHE_FN):
              try:
                  angleInfo = None
                  cache = shelve.open(self.CACHE_FN)
                  self.num_res, sim_chis = cache[key]
                  cache.close()
                  found_in_cache = True
                  print "Loaded '%s' from the cache" % key, flush()
              except Exception,e:
                  print e
                  pass # if the data isn't found in the cache, then lod it

      # if it wasn't found in the cache, load it
      if not found_in_cache:
          print "Loading trajectory '%s'" % pdb_traj.traj_fn, flush()

          pdb_num = 0
          for pdb in pdb_traj.get_next_pdb():
              if pdb_num >= self.num_structs: continue

              pdb_chis = pdb.calc_chis()
              pdb_phipsiomega = pdb.calc_phi_psi_omega()
              if (pdb_num+1) % 100 == 0:
                  print "Loaded pdb #%d" % (pdb_num), flush()

              if sim_chis == None:
                  self.num_res = pdb_chis.shape[0]
                  sim_chis = zeros((self.num_structs, self.num_res,6), float64)
              sim_chis[pdb_num, :, 0:2] = pdb_phipsiomega[:,0:2]
              sim_chis[pdb_num, :, 2:] = pdb_chis
              #print " shape phipsiomega:"+str(shape(pdb_phipsiomega))+" shape chis: "+str(pdb_chis.shape[0])
              pdb_num += 1
          #print "sim chis", sim_chis

      # Store the angle data into the cache
      if cache_to_disk and not found_in_cache:
          cache = shelve.open(self.CACHE_FN)
          cache[key] = [self.num_res, sim_chis]
          cache.close()

      if self.all_chis == None:
          self.all_chis = zeros((self.num_sims, self.num_structs, self.num_res,6), float64)
      self.all_chis[sequential_sim_num, :, : , :] = sim_chis[:, :, :]

   def get_angles(self, sequential_sim_num, sequential_res_num, res_name, backbone_only, phipsi):
       chi_nums = []
       if phipsi == 2: chi_nums += [0,1] # which angles to use
       if backbone_only == 0: chi_nums += range(2,NumChis[res_name]+2)
       #print chi_nums
       #print self.all_chis

       curr_angles = zeros((len(chi_nums), self.all_chis.shape[1]))
       for sequential_chi_num, chi_num in zip(range(len(chi_nums)), chi_nums):
           curr_angles[sequential_chi_num, :] = self.all_chis[sequential_sim_num, :, int(sequential_res_num)-1, chi_num]
       return curr_angles, curr_angles.shape[1]

#========================================
#FUNCTION DEFINITIONS
#========================================
#Function definition to read an xvg file
def readxvg(filename,skip):
   """Read desired xvg file; strip headers and return data as array. First column of array is times of data points; remaining columns are the data. Should properly truncate the end of the data file if any of the lines are incomplete.
INPUT: Name or path of xvg file to read.
RETURN: As a tuple:
(1) An LxN data array containing the data from the file, less the header and any aberrant lines from the end (aberrant in the sense of truncated or not following the pattern of the rest of the lines). N is the number of columns in the xvg file, and L the number of lines. It is up to the user to interpret these.
Note that units are as in xvg file (normall kJ/mol for energies from GROMACS)
(2) The title as read from the xvg file
"""

   #Read input data   
   fil=open(filename,'r');
   inlines=fil.readlines()
   fil.close()

   #Slice off headers
   #Find header lines beginning with @ or #.
   headerline=re.compile(r'[@#].*')
   match=True
   linenum=0
   title=''
   while (match):
     m=headerline.match(inlines[linenum])
     if not m:
        match=False
     else:
        #obtain title
        if inlines[linenum].find('title')>-1:
           tmp=inlines[linenum].split() 
           title=tmp[2]+' '+tmp[3]
        #Go to next line
        linenum+=1
   #slice off headers
   inlines=inlines[linenum:]

   #Detect how many fields on each line in body of xvg file. 
   numfields=len(inlines[0].split())

   #Length (including any aberrant lines at the end) 
   inlength=len(inlines)

   #Array to store data
   dataarray=zeros((int(inlength/skip),numfields),float64)

   skiplines=0
   #Read data into array
   for i in range(int(inlength/skip)):
      entries=inlines[i*skip].split()
      #Make sure find expected number of entries on line
      tmpentries=len(entries)
      if tmpentries!=numfields:
        print "Found %(tmpentries)s on line %(i)s; expected %(numfields)s. Skipping line and continuing." % vars()
        skiplines+=1
      elif entries[1]=='nan':
        #Do a bit of checking also for corrupted data as in the case of corrupted trajectories
        #which sometimes give nan on this step.
        skiplines+=1
        print "Found some 'nan' entries on line %(i)s. Skipping." % vars()
      else:
        #Store data to data array, in packed format
        for j in range(numfields):
           dataarray[i-skiplines][j]=float(entries[j])

   #Last (skiplines) of dataarray will be empty, so pack data array
   dataarray=resize(dataarray,(int(inlength/skip)-skiplines,numfields))
 
   return (dataarray,title)

# output the diagonal elements of a matrix
def output_diag(myfilename,mymatrix,rownames):
   #outputs only diagonal
   myfile = open(myfilename,'w')
   for row_num, row_name in zip(range(len(rownames)), rownames):
      myfile.write(row_name + " ")
      for col_num, col_name in zip(range(len(rownames)), rownames):
         if(row_name == col_name):
            myfile.write(str(mymatrix[row_num,col_num]))
            myfile.write("\n")
   myfile.close()

def output_entropy(myfilename,mylist):
   myfile = open(myfilename,'w')
   for i in range(mylist.shape[0]):
       myfile.write(str(mylist[i]))
       myfile.write("\n")
   myfile.close()

def output_value(myfilename,myvalue):
   myfile = open(myfilename,'w')
   myfile.write(str(myvalue))
   myfile.write("\n")
   myfile.close()

# output the elements of a matrix in string formatting, optionally zeroing the diagonal terms
def output_matrix(myfilename,mymatrix,rownames,colnames, zero_diag=False):
   myfile = open(myfilename,'w')

   for col_num, col_name in zip(range(len(colnames)), colnames):
      myfile.write(col_name + " ")
   myfile.write("\n")
   for row_num, row_name in zip(range(len(rownames)), rownames):
      myfile.write(row_name + " ")
      for col_num, col_name in zip(range(len(colnames)), colnames):
         if col_num == row_num and zero_diag:
            myfile.write(str(0))
         else:
            #print row_num, col_num, mymatrix[row_num,col_num]
            myfile.write(str(mymatrix[row_num,col_num]))
         myfile.write(" ")
      myfile.write("\n")
   myfile.close()

def output_matrix_chis(myfilename,mymatrix,rownames,colnames, nchi=6, zero_diag=False):
   myfile = open(myfilename,'w')
   #print "shape of matrix to ouput:"+str(mymatrix.shape)+"\n"
   for col_num, col_name in zip(range(len(colnames)), colnames):
     for col_chi in range(nchi):
      myfile.write(col_name + "_" +str(col_chi) + " ")
   myfile.write("\n")
   for row_num, row_name in zip(range(len(rownames)), rownames):
     for row_chi in range(nchi):
      myfile.write(row_name + "_" + str(row_chi) + " ")
      for col_num, col_name in zip(range(len(colnames)), colnames):
        for col_chi in range(nchi):  
         if col_num == row_num and row_chi == col_chi and zero_diag:
            myfile.write(str(0))
         else:
            #print row_num, col_num, mymatrix[row_num,col_num]
            myfile.write(str(mymatrix[row_num,col_num,row_chi,col_chi]))
         myfile.write(" ")
      myfile.write("\n")
   myfile.close()


def read_matrix_chis(myfilename, nchi=6, zero_diag=False):
   rownames = []
   colnames = []
   myfile = open(myfilename,'r')
   inlines = myfile.readlines()
   myfile.close()
   reschis = inlines[0].split()
   mymatrix = zeros((int(len(inlines[1:]) / nchi), int((len(reschis))/nchi),6,6),float64)
   #print mymatrix.shape
   for myname_num in reschis:
       (thisname, thisnum) = myname_num.split('_')
       if int(thisnum) == 0:
           colnames.append(thisname)
   #print colnames
   #print len(colnames)
   for row_num in range(int(len(inlines[1:]))):
       thisline = inlines[row_num + 1]
       thislinedata = thisline.split()
       (thisname, row_chi) = thislinedata[0].split('_')
       res_num = int(floor(row_num / nchi))
       row_chi = int(row_chi) #convert string value to integer
       thislinenums = map(float, thislinedata[1:]) #does this need to be float64 or another double precision thingy?
       #print thislinenums
       thislinearray = array(thislinenums,float64)
       #print thislinearray.shape
       if row_chi == 0:
           rownames.append(thisname)
       for col_num in range(len(colnames)):
           for col_chi in range(nchi):
               #print "name: "+str(thisname)+" chi: "+str(row_chi)+ " row_num: "+str(row_num)+" row_chi: "+str(row_chi)+ " col_num: "+str(col_num)+" col_chi: "+str(col_chi)+"\n"
               mymatrix[res_num,col_num,row_chi,col_chi] = float64(thislinearray[col_num*nchi + col_chi])
   #print rownames
   return mymatrix, rownames, colnames


def read_res_matrix(myfilename):
   rownames = []
   colnames = []
   myfile = open(myfilename,'r')
   inlines = myfile.readlines()
   myfile.close()
   res = inlines[0].split()
   mymatrix = zeros((int(len(inlines[1:])), int(len(res))),float64)
   #print mymatrix.shape
   for myname_num in res:
       colnames.append(thisname)
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




def compute_CA_dist_matrix(reslist,pdb_obj):
    CA_dist_matrix = zeros((len(reslist),len(reslist)),float32)
    
    for res_ind1, myres1 in zip(range(len(reslist)), reslist):
        for res_ind2, myres2 in zip(range(res_ind1, len(reslist)), reslist[res_ind1:]):
            print "\n#### Working on residues %s and %s (%s and %s):" % (myres1.num, myres2.num, myres1.name, myres2.name) 
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

usage="%prog [-s mypdb.pdb] resfile  # where resfile is in the format <1-based-index> <aa type> <res num>"
parser=OptionParser(usage)
parser.add_option("-s", "--structure", default=None, type="string", help="filenames to load PDB from")
parser.add_option("-m", "--matrix",    default=None, type="string", help="res x res mutinf matrix to read from")
(options,args)=parser.parse_args()

if len(filter(lambda x: x==None, (options.structure))) == 1:
    parser.error("ERROR exactly one of --structure must be specified")

resfile_fn = args[0]



run_params = RunParameters(num_structs=1,resfile_fn=resfile_fn)
pdb_lines = commands.getoutput("egrep '^ATOM ' %s" % options.structure).split("\n")
pdb_obj = PDBlite.PDB(pdb_lines,fn=options.structure)
print pdb_obj
reslist = load_resfile(run_params,load_angles=False)
#print reslist
CA_dist_matrix = compute_CA_dist_matrix(reslist,pdb_obj)
name_num_list=[]
for res in reslist: name_num_list.append(res.name + str(res.num))
output_matrix(resfile_fn+"_CA_dist_matrix.txt",CA_dist_matrix,name_num_list,name_num_list)

if(options.matrix =! None):
    mutinf_res_matrix = read_res_matrix(options.matrix)
    # now, filter res matrix at various fractions of kT
    # start with (S/k * kT) >= 1/2 kT
    half_kT_mutinfs = mutinf_res_matrix.copy()
    half_kT_mutinfs(half_kT_mutinfs < 0.5 * (8.314*300/(4.184*1000)) ) = 0
    #next, take (S/k * kT) >= 1/4 kT
    quarter_kT_mutinfs = mutinf_res_matrix.copy()
    quarter_kT_mutinfs(quarter_kT_mutinfs < 0.25 * (8.314*300/(4.184*1000)) ) = 0
    #next, take (S/k * kT) >= 0.1 kT
    tenth_kT_mutinfs = mutinf_res_matrix.copy()
    tenth_kT_mutinfs(tenth_kT_mutinfs < 0.25 * (8.314*300/(4.184*1000)) ) = 0

    output_matrix(resfile_fn+"_tenth_kT_mutinfs.txt",tenth_kT_mutinfs,name_num_list,name_num_list)
    #now, grab distances for values that pass

    half_kT_CA_dist_matrix = CA_dist_matrix.copy()
    half_kT_CA_dist_matrix(half_kT_mutinfs == 0) = 0
    output_matrix(resfile_fn+"_half_kT_CA_dist_matrix.txt",half_kT_mutinfs,name_num_list,name_num_list)

    quarter_kT_CA_dist_matrix = CA_dist_matrix.copy()
    quarter_kT_CA_dist_matrix(quarter_kT_mutinfs == 0) = 0
    output_matrix(resfile_fn+"_quarter_kT_CA_dist_matrix.txt",quarter_kT_mutinfs,name_num_list,name_num_list)

    tenth_kT_CA_dist_matrix = CA_dist_matrix.copy()
    tenth_kT_CA_dist_matrix(tenth_kT_mutinfs == 0) = 0
    output_matrix(resfile_fn+"_tenth_kT_CA_dist_matrix.txt",tenth_kT_mutinfs,name_num_list,name_num_list)




#END

