#!/usr/bin/python
#Python script to analyze dihedrals. Usage: "python dihedral_analyze.py (list of files to plot). D. Mobley, 5/16/07"
from numpy import zeros, arange, outer, inner, float32, int16, sqrt, floor, mean, std, log, add, average, nonzero, cov, transpose, resize, array, reshape, max, array2string, sum, histogram2d, logical_and, random, sort, searchsorted, swapaxes, shape
import re, os, sys, os.path, time, shelve
from optparse import OptionParser
from scipy import weave
from scipy import stats as stats
from scipy.weave import converters
import PDBlite

SMALL = 0.000000001
CACHE_TO_DISK = True
VERBOSE = 1
OFF_DIAG = 1 #set to "1" to look at off-diagonal terms
OUTPUT_DIAG = 0  #set to "1" to output flat files for diagonal matrix elements
MULTIPLE_1D_BINS = 2 #set to "2" for double the number of 1-D bins for 1-D entropy only
SPEED = "vfast" # set to "vfast" to use the vfast_cov(); "fast" to use fast_cov(); "*" to use cov()
print "Running at speed: ", SPEED

# from: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc
                
# in seconds
class Timer:
    def __init__(self): self._start = time.time()
    def elapsed(self):  return time.time() - self._start
    def __str__(self):  return "%0d:%02d" % (self.elapsed()/60, self.elapsed()%60)

# convert an array to a pretty string
def arr2str(a): return array2string(a, max_line_width=160, precision=1)

# calculates the covariance matrix over the rows
cov_mat = None
def fast_cov(m):
   global cov_mat

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float32)
   X -= X.mean(axis=0)

   if cov_mat == None: cov_mat = zeros((ny, ny), float32)
   cov_mat[:,:] = 0

   code = """
      // weave
      float inv_nx, cov;

      inv_nx = 1./float(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         for (int y2=y1; y2<ny; ++y2) {
            cov=0;
            for (int x=0; x<nx; ++x) {
               cov += X(x,y1) * X(x,y2);
            }
            cov_mat(y1,y2) = cov*inv_nx;
            cov_mat(y2,y1) = cov_mat(y1,y2);
         }
      }
      //printf("fast_cov max=%f\\n", max);
   """
   weave.inline(code, ['X', 'cov_mat', 'nx', 'ny'],
                type_converters = converters.blitz, compiler = 'gcc')
   return cov_mat


# calculates the covariance matrix over the rows
def vfast_cov(m):
   global cov_mat

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float32)
   X -= X.mean(axis=0)

   if cov_mat == None: cov_mat = zeros((ny, ny), float32)
   cov_mat[:,:] = 0
   
   code = """
      // weave
      float inv_nx, cov;
      float *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./float(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         pCovy1 = cov_mat+y1*ny;
         pXy1 = X+y1;
         for (int y2=y1; y2<ny; ++y2) {
            cov=0;
            pXy2 = X+y2;
            for (int x=0; x<nx; ++x) {
               //cov += X(x,y1) * X(x,y2);
               //cov += *(X+x*ny+y1) * *(X+x*ny+y2);               
               tmp=x*ny;
               cov += *(pXy1+tmp) * *(pXy2+tmp);
            }
            //cov_mat(y1,y2) = cov*inv_nx;
            //cov_mat(y2,y1) = cov_mat(y1,y2);
            *(pCovy1+y2) = cov*inv_nx;
            *(cov_mat+y2*ny+y1) = *(pCovy1+y2);
         }
      }
   """
   weave.inline(code, ['X', 'cov_mat', 'nx', 'ny'],
                #type_converters = converters.blitz,
                compiler = 'gcc')
   return cov_mat

#cov_mat_for1D = None
def vfast_cov_for1D(m):
   global cov_mat_for1D

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float32)
   X -= X.mean(axis=0)

   cov_mat_for1D = zeros((ny, ny), float32)
   cov_mat_for1D[:,:] = 0
   
   code = """
      // weave
      float inv_nx, cov;
      float *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./float(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         pCovy1 = cov_mat_for1D+y1*ny;
         pXy1 = X+y1;
         for (int y2=y1; y2<ny; ++y2) {
            cov=0;
            pXy2 = X+y2;
            for (int x=0; x<nx; ++x) {
               //cov += X(x,y1) * X(x,y2);
               //cov += *(X+x*ny+y1) * *(X+x*ny+y2);               
               tmp=x*ny;
               cov += *(pXy1+tmp) * *(pXy2+tmp);
            }
            //cov_mat_for1D(y1,y2) = cov*inv_nx;
            //cov_mat_for1D(y2,y1) = cov_mat_for1D(y1,y2);
            *(pCovy1+y2) = cov*inv_nx;
            *(cov_mat_for1D+y2*ny+y1) = *(pCovy1+y2);
         }
      }
   """
   weave.inline(code, ['X', 'cov_mat_for1D', 'nx', 'ny'],
                #type_converters = converters.blitz,
                compiler = 'gcc')
   return cov_mat_for1D

cov_mat_for1D_ind = None
def vfast_cov_for1D_ind(m):
   global cov_mat_for1D_ind

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float32)
   X -= X.mean(axis=0)

   if cov_mat_for1D_ind == None: cov_mat_for1D_ind = zeros((ny, ny), float32)
   cov_mat_for1D_ind[:,:] = 0
   
   code = """
      // weave
      float inv_nx, cov;
      float *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./float(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         pCovy1 = cov_mat_for1D_ind+y1*ny;
         pXy1 = X+y1;
         for (int y2=y1; y2<ny; ++y2) {
            cov=0;
            pXy2 = X+y2;
            for (int x=0; x<nx; ++x) {
               //cov += X(x,y1) * X(x,y2);
               //cov += *(X+x*ny+y1) * *(X+x*ny+y2);               
               tmp=x*ny;
               cov += *(pXy1+tmp) * *(pXy2+tmp);
            }
            //cov_mat_for1D_ind(y1,y2) = cov*inv_nx;
            //cov_mat_for1D_ind(y2,y1) = cov_mat_for1D_ind(y1,y2);
            *(pCovy1+y2) = cov*inv_nx;
            *(cov_mat_for1D_ind+y2*ny+y1) = *(pCovy1+y2);
         }
      }
   """
   weave.inline(code, ['X', 'cov_mat_for1D_ind', 'nx', 'ny'],
                #type_converters = converters.blitz,
                compiler = 'gcc')
   return cov_mat_for1D_ind


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
        return "Num sims = %d; Runs being used = %s; Binwidth = %.1f; Permutations = %d; Sigma = %.1f" % \
               (d['num_sims'], " ".join(map(str, d['which_runs'])), d['binwidth'], d['permutations'], d['sigma'])

    def get_logfile_prefix(self):
        return "%s-nsims%d-structs%d-bin%d" % (os.path.basename(self.resfile_fn), self.num_sims, self.num_structs, int(self.binwidth))

# Class to hold all the dihedral angle info for a set of trajectories
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
              except Exception,e:
                  print e
                  pass # if the data isn't found in the cache, then lod it

      # if it wasn't found in the cache, load it
      if not found_in_cache:
          print "Loading trajectory '%s'" % pdb_traj.traj_fn
          sys.stdout.flush()

          pdb_num = 0
          for pdb in pdb_traj.get_next_pdb():
              if pdb_num >= self.num_structs: continue

              pdb_chis = pdb.calc_chis()
              pdb_phipsiomega = pdb.calc_phi_psi_omega()
              if (pdb_num+1) % 100 == 0:
                  print "Loaded pdb #%d" % (pdb_num)
                  sys.stdout.flush()

              if sim_chis == None:
                  self.num_res = pdb_chis.shape[0]
                  sim_chis = zeros((self.num_structs, self.num_res, 6), float32)
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
          self.all_chis = zeros((self.num_sims, self.num_structs, self.num_res, 6), float32)
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
def readxvg(filename):
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
   dataarray=zeros((inlength,numfields),float32)

   skiplines=0
   #Read data into array
   for i in range(inlength):
      entries=inlines[i].split()
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
   dataarray=resize(dataarray,(inlength-skiplines,numfields))
 
   return (dataarray,title)


def bintouple(angle1,angle2,binwidth):
   bin1 = int(floor((angle1-0.00001 + 180) / binwidth))
   bin2 = int(floor((angle2-0.00001 + 180) / binwidth))
   return [bin1, bin2]

def binsingle(angle,inv_binwidth):
   if angle == -180: angle = 180
   return int(floor((angle-0.00001 + 180)*inv_binwidth)) #so we don't get an overshoot if angle is exactly 180

def binsingle_adaptive(angle,inv_binwidth):
    #print "rank: "+str(angle)+" binwidth: "+str(1.0/inv_binwidth)+" bin: "+str(int(floor(angle*inv_binwidth)))
    return int(floor(angle*inv_binwidth)) #here "angle" is a rank-order for the angle over sum(numangles)
   
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

# counts has shape [num_sims x num_angles]
def calc_entropy(counts, num_sims, nchi, bootstrap_choose=0, calc_variance=True, entropy = None, var_ent = None):
    av_counts=average(counts,axis=0)
    #print av_counts
    #assert(entropy != None and var_ent != None)
    log_av_counts = log(av_counts + SMALL)
    #print -sum(av_counts * ( log_av_counts),axis=1)
    #print entropy
    entropy[:] = -sum(av_counts * ( log_av_counts),axis=1) #this is not needed for off-diagonal terms as these marginal entropies cancel out with terms in the mutual information correction
    if(calc_variance):
        deriv_vector = - (1 + log_av_counts)
    #print "chi: "+str(mychi)+" bin: "+str(b)+" av_counts: " +str(av_counts)+"\n"+" entropy: "+str(-av_counts * ( log_av_counts))+" unc_av: "+str(unc_counts)+" var_ent this chi: "+str((- (1 + log_av_counts)) * (- (1 + log_av_counts)) * (unc_counts * unc_counts))
    if(num_sims > 1 and calc_variance):
        for mychi in range(nchi):   
            bins_cov = vfast_cov_for1D(counts[mychi,:,:]) #observations are in the rows, conditions(variables) are in the columns
            #print bins_cov.shape, deriv_vector.shape
            #now we will gather all terms, diagonal and cross-terms, in the variance of joint entropy using matrix multiplication instead of explicit loop
            # var = A * Vij * T(A), where A={dy/dxj}|x=mu
            if VERBOSE >= 2:
                print av_counts.shape, deriv_vector.shape, bins_cov.shape, transpose(deriv_vector).shape
            var_ent[mychi] = inner(deriv_vector,inner(bins_cov,transpose(deriv_vector)))
         
    return entropy, var_ent

# Calculate the MI between two angle distributions, x & y.
# First calculate the vector (over all bins) of means of the joint probabilities, Pij.
# Then calculate the vectors (over all bins) of Pi*Pj and Pi+Pj.
# Now MI(x,y) = sum over all bins [ Pij * log(Pij/(Pi*Pj)) ].
# The variance of MI = d*C*d'    where d={dy/dxj}|x=mu = 1 + log(Pij/(Pi*Pj)) - (Pi+Pj) * Pij/(Pi*Pj)
#     and C is the matrix of the covariances Cij = sum over angles i [ (xi-mu(x)) * (yi-mu(y)) ]/(N-1)
pop_matrix = U = logU = None
def calc_mutinf(chi_pop_hist1, chi_pop_hist2, bins1, bins2, num_sims, nbins, numangles, calc_variance=True,bootstrap_choose=0,permutations=0):
    global pop_matrix, U, logU

    # allocate the matrices the first time only for speed
    if pop_matrix == None:
       pop_matrix = zeros((permutations + 1 , num_sims, nbins*nbins), float32)
       U = zeros((permutations + 1, nbins*nbins), float32)
       logU = zeros((permutations + 1, nbins*nbins),float32)   
    pop_matrix[:,:,:] = 0
    
    code = """
      // weave
      int mynumangles = 0;
      for (int permut=0; permut < permutations + 1; permut++) {
       for (int simnum=0; simnum < num_sims; simnum++) {
          mynumangles = *(numangles + simnum);
          for (int anglenum=0; anglenum< mynumangles; anglenum++) {            
             *(pop_matrix + permut*nbins*nbins*num_sims + nbins*nbins*simnum
              +   (*(bins1 + simnum*mynumangles + anglenum))*nbins
              +   (*(bins2 + permut*num_sims*mynumangles + simnum*mynumangles + anglenum))) += 1;
          }
        }
       }
      """
    weave.inline(code, ['num_sims', 'numangles', 'nbins', 'bins1', 'bins2', 'pop_matrix','permutations'],
                 #type_converters = converters.blitz,
                 compiler = 'gcc')


    #code = """
    #  // weave
    #  int mynumangles = 0;
    #  for (int simnum=0; simnum < num_sims; simnum++) {
    #     mynumangles = *(numangles + simnum);
    #     for (int anglenum=0; anglenum< mynumangles; anglenum++) {            
    #        *(pop_matrix + nbins*nbins*simnum
    #         +   (*(bins1 + simnum*mynumangles + anglenum))*nbins
    #         +   (*(bins2 + simnum*mynumangles + anglenum))) += 1;
    #     }
    #  }
    #  """
    #weave.inline(code, ['num_sims', 'numangles', 'nbins', 'bins1', 'bins2', 'pop_matrix'],
    #             #type_converters = converters.blitz,
    #             compiler = 'gcc')
    
    #if VERBOSE >=2:
    #print "pop_matrix  before normalization:"
    #print str(sum(pop_matrix[0,:,:],axis=0))
    
    for simnum in range(num_sims): pop_matrix[:,simnum,:] /= numangles[simnum]
    if bootstrap_choose > 1:
        print "map pop_matrix from nsims into nsims choose simnum sims"
    
    ##################################
    #output pop matrix for debugging
    #mybinlist = []
    #for mybin in range(nbins):
    #   mybinlist.append(str(mybin))
    #mymatrix = reshape(pop_matrix,(num_sims,nbins,nbins))
    #temp juryrig -- not robust for runs of different sim length
    #myfilename = self.name + str(self.num) + str(mychi1+1) + str(mychi2+1)+"pophist.txt"
    #output_matrix(myfilename,sum(mymatrix * self.numangles[0],axis=0),mybinlist,mybinlist)
    #done with output pop matrix
    #print pop_matrix[0,:]
    #may want to make a plot of this
    ###################################

    #now calculate mutual information by integrating over bins
    #do calcuations in 1-D: accumulate entropy, convert deriv_matrix back to 1-D
    #print chi_pop_hist1 * sum(numangles)
    #print chi_pop_hist2 * sum(numangles)
    #print sum(outer(chi_pop_hist1,chi_pop_hist2)[:,:],axis=1)
    pxipyj_flat = outer(chi_pop_hist1,chi_pop_hist2).flatten()
    pxipyj_flat = resize(pxipyj_flat,(permutations + 1,(pxipyj_flat.shape)[0]))
    av_joint_pop = average(pop_matrix,axis=1)
    #print "av_joint_pop\n"
    #print av_joint_pop[0,:]
    #print "pxipjy_flat\n"
    #print (pxipyj_flat[:,:] + SMALL)
    U[:,:] = 0
    logU[:,:] = 0
    U = av_joint_pop / (pxipyj_flat + SMALL)
    #print "U"
    #print U[0,:]
    logU = log(U + SMALL)
    #print "logU"
    #print logU[0,:]
    #for mypermut in range(1):
        #nonzero_bins = logical_and((av_joint_pop[mypermut,:] != 0), (pxipyj_flat != 0))

        
        #U[mypermut, nonzero_bins] = av_joint_pop[mypermut, nonzero_bins] / pxipyj_flat[nonzero_bins]
        #logU[mypermut, nonzero_bins] = log(U[mypermut, nonzero_bins])
    #print (av_joint_pop * logU)[0,:]
    #print logU.shape

    mutinf_thisdof = sum(av_joint_pop * logU,axis=1)

    
    #av_joint_pop = average(pop_matrix,axis=1)
    #print av_joint_pop
    #pxipyj_flat = outer(chi_pop_hist1,chi_pop_hist2).flatten() + SMALL
    #resize(pxipyj_flat,(permutations + 1,(pxipyj_flat.shape)[0]))
    #print pxipyj_flat.shape
    #print pxipyj_flat
    #nonzero_bins = logical_and((av_joint_pop != 0), (pxipyj_flat != 0))

    #U[:,:] = 0
    #logU[:,:] = 0
    #U = av_joint_pop / pxipyj_flat + SMALL
    #logU = log(U)
    #mutinf_thisdof = sum(av_joint_pop * logU,axis=1)

    #print mutinf_thisdof
    #myfilename = self.name + str(self.num) + str(mychi1+1) + str(mychi2+1)+"mutentbeforesum.txt"
    #mymatrix = reshape((av_joint_pop * logU),(nbins,nbins))
    #output_matrix(myfilename,mymatrix,mybinlist,mybinlist)

    if calc_variance == False:
       #if VERBOSE: print "   mutinf = %.5f," % (mutinf_thisdof),
       return mutinf_thisdof, 0

    # print out the number of nonzero bins
    if VERBOSE >=2:
     #  num_nonzero_jointpop, num_nonzero_pxipyj = len(nonzero(av_joint_pop)[0]), len(nonzero(pxipyj_flat)[0]), 
     #  print "   nonzero bins (tot=%d): jointpop=%d, pxipyj=%d, combined=%d" % (nbins*nbins, num_nonzero_jointpop, num_nonzero_pxipyj, len(nonzero_bins[nonzero_bins==True]))
       print "   mutinf this degree of freedom "+str(mutinf_thisdof)
    # calculate the variance of the mutual information
    pxi_plus_pyj_flat = add.outer(chi_pop_hist1, chi_pop_hist2).flatten()
    deriv_vector = 1 + logU[0,:] - (pxi_plus_pyj_flat[0,:]) * U
    if(num_sims > 1):
       if SPEED == "vfast":
          bins_cov = vfast_cov(pop_matrix[0,:,:]) #observations are in the rows, conditions(variables) are in the columns
       elif SPEED == "fast":
          bins_cov = fast_cov(pop_matrix[0,:,:]) #observations are in the rows, conditions(variables) are in the columns
       else:
          bins_cov = cov(pop_matrix[0,:,:],rowvar=0) #observations are in the rows, conditions(variables) are in the columns
          if VERBOSE > 2:
              print bins_cov
       #print deriv_vector
       #print bins_cov
       #now we will gather all terms, diagonal and cross-terms, in the variance of joint entropy using matrix multiplication instead of explicit loops
       # var = A * Vij * T(A), where A={dy/dxj}|x=mu
       var_mi_thisdof = inner(deriv_vector,inner(bins_cov,transpose(deriv_vector)))
    else:
       var_mi_thisdof = 0

    return mutinf_thisdof, var_mi_thisdof

### calculation of independent information using permutations
independent_mutinf_thisdof = None
def calc_independent_mutinf(chi_pop_hist1, chi_pop_hist2, bins1, bins2, num_sims, nbins, numangles, permutations, bootstrap_choose):
   if permutations == 0: return (0, 0)
   global independent_mutinf_thisdof 

   # permutations of the distribution are used for independent mutual information, which when subtracted from mutual information, gives the excess mutual information
   if(independent_mutinf_thisdof == None):
       independent_mutinf_thisdof = zeros((permutations,1),float32)
   for i in range(permutations):
      independent_mutinf_thisdof[i,0] = calc_mutinf(chi_pop_hist1, chi_pop_hist2, bins1, bins2scrambled, num_sims, nbins, numangles, calc_variance=False,bootstrap_choose=bootstrap_choose)[0]

   #if VERBOSE >= 2: print "   ind mutinfs = %s" % " ".join(map(lambda x:"%.3f"%x, independent_mutinf_thisdof))
   avg_independent_mutinf_thisdof = sum(independent_mutinf_thisdof)/permutations
   var_independent_mutinf_thisdof = (vfast_cov_for1D_ind(independent_mutinf_thisdof))[0,0]
   return avg_independent_mutinf_thisdof, var_independent_mutinf_thisdof

def calc_excess_mutinf(chi_pop_hist1, chi_pop_hist2, bins1, bins2, num_sims, nbins, numangles, sigma, permutations, bootstrap_choose, calc_variance=False):
    if VERBOSE >= 2:
       print "   chi_pop_hist1", int16(chi_pop_hist1 * sum(numangles)), "  length:", len(chi_pop_hist1)
       print "   chi_pop_hist2", int16(chi_pop_hist2 * sum(numangles)), "  length:", len(chi_pop_hist2)

    mutinf_tot_thisdof, var_mi_thisdof = calc_mutinf(chi_pop_hist1,chi_pop_hist2, bins1, bins2, num_sims, nbins, numangles, calc_variance=False, bootstrap_choose=bootstrap_choose, permutations=permutations) 
    #independent_mutinf_thisdof, var_ind_mi_thisdof = calc_independent_mutinf(chi_pop_hist1,chi_pop_hist2, bins1, bins2, num_sims, nbins, numangles, permutations, bootstrap_choose) 
    #mutinf_tot_thisdof = mutinf_thisdof.copy()
    var_mi_thisdof_old = var_mi_thisdof
    independent_mutinf_thisdof = average(mutinf_tot_thisdof[1:])
    excess_mutinf_thisdof = mutinf_tot_thisdof[0] - independent_mutinf_thisdof
    test_stat = 0
    mycutoff = 0
    sigtext = " "
    #if mutinf_thisdof[0] < 0:
    #    mutinf_thisdof[0] = var_mi_thisdof = test_stat = 0
    #    sigtext = " *not significant*"
    #    mycutoff = 1
    #    print "   mutinf/ind_mutinf = %.3f %.3f (sd= %.3f ; S= %.3f ; cut= %.3f )   %s" % (mutinf_thisdof, independent_mutinf_thisdof, sqrt(var_mi_thisdof_old), test_stat, mycutoff, sigtext)
    #return excess_mutinf_thisdof, var_mi_thisdof, test_stat
    #else:
    # test_stat = mycutoff = 0
    # sigtext = " # "
    # if(permutations > 0):
    #    ind_gamma_scale = independent_mutinf_thisdof * independent_mutinf_thisdof / var_ind_mi_thisdof
    #    ind_gamma_shape = var_ind_mi_thisdof / independent_mutinf_thisdof
    #    if(num_sims == 1):
    #        mydist = stats.gamma
    #        mycutoff = mydist.ppf(0.95,ind_gamma_shape,loc=0,scale=ind_gamma_scale)
    #        test_stat = mutinf_tot_thisdof
    #        #S = stats.gamma.sf(mutinf_tot_thisdof,ind_gamma_scale,loc=0,scale=ind_gamma_shape)        
    #    else:
    #        excess_gamma_scale = mutinf_tot_thisdof * mutinf_tot_thisdof / var_mi_thisdof
    #        excess_gamma_shape = var_mi_thisdof / mutinf_tot_thisdof
    #        mydist=stats.f
    #        test_stat = mutinf_tot_thisdof / independent_mutinf_thisdof
    #        mycutoff = mydist.ppf(0.95, 2 * num_sims * excess_gamma_scale, 2 * permutations * ind_gamma_scale)
    #        #S = mycutoff * independent_mutinf_thisdof

    
    # if (permutations > 0 and test_stat <= mycutoff):
    #		mutinf_thisdof = var_mi_thisdof = 0
    #            sigtext = " *not significant*"
    
    print "   mutinf/ind_mutinf = %.3f %.3f (sd= %.3f ; S= %.3f ; cut= %.3f )   %s" % (mutinf_tot_thisdof[0], independent_mutinf_thisdof, 0, 0, 0, sigtext)
    
    #use sigma as a cutoff for 'significant' mutents if we didn't do any permutations
    #print "   mutinf/ind_mutinf = %.3f %.3f (sd= %.3f ; S= %.3f )" % (mutinf_thisdof, independent_mutinf_thisdof, sqrt(var_mi_thisdof), S),
    #if S>2.132: print "#",
    #try 1.553 for 90%, 2.132 for 95% for 1-sided Student's t-test with n-1 = 4 degrees of freedom, S = (MI - mu) / ( s / root n), where s is sample stdev
    #if (var_mi_thisdof * sigma >= mutinf_thisdof * mutinf_thisdof):
    #if ( S < 2.132):
     #if (S < mycutoff):
     #  mutinf_thisdof = var_mi_thisdof = 0
     #  print " *not significant*",
     #  print

    return excess_mutinf_thisdof, var_mi_thisdof, test_stat

# Number of chi angles per residue
NumChis = { "ALA":0, "CYS":1, "CYN":1, "CYX":1, "ASP":2, "GLU":3, "GLH":3, "PHE":2, "GLY":0, "HIS":2,
            "HIP":2, "HIE":2, "HID":2, "ILE":2, "LYS":4, "LYP":4, "LEU":2, "MET":3, "ASN":2, "GLN":3,
            "PRO":1, #GF: can't find def for chi2; rosetta uses chi1 only,
            "ARG":4, "SER":1, "THR":1, "VAL":1,"TRP":2, "TYR":2, "CTH":2}


class ResidueChis:
   name = 'XXX'
   num = 0
   nchi = 0
   angles = []
   rank_order_angles = []
   angles_complex = [] #for spectral density calculations
   chi_pop_hist = []
   chi_var_pop = []
   bins = []
   entropy = 0
   var_ent = 0
   inv_numangles = []
   counts = []
   # load angle info from xvg files
   def _load_xvg_data(self, basedir, num_sims, max_angles, chi_dir = "/dihedrals/g_chi/"):
      myname = self.name
      mynumchis = NumChis[myname]
      shifted_angles = zeros((self.nchi,num_sims,max_angles),float32)
      #shifted_angles[:,:,:] = -999 #a null value other than zero
      #weird fix for residue type "CYS2" in Gromacs
      if myname == "CYS": myname += "2"
      assert num_sims == len(self.which_runs)
          
      for chi_num in range(self.nchi):
         #print "Chi:"+str(mychi+1)+"\n"
         for sequential_sim_num in range(num_sims):
            sim_index_str = str(self.which_runs[sequential_sim_num])
            if(chi_num < mynumchis):
               xvg_fn = basedir+"run"+sim_index_str+chi_dir+"chi"+str(chi_num+1)+myname+self.num+".xvg"
            if(chi_num == mynumchis):
               xvg_fn = basedir+"run"+sim_index_str+chi_dir+"phi"+myname+self.num+".xvg"
            if(chi_num == mynumchis + 1):
               xvg_fn = basedir+"run"+sim_index_str+chi_dir+"psi"+myname+self.num+".xvg"
               
            # gf: gromacs is changing some of my residue names when I use it to read in a PDB trajectory
            if not os.path.exists(xvg_fn):
               if myname == "LYS":
                  for myname_alt in ("LYSH",):
                     xvg_fn = basedir+"run"+sim_index_str+chi_dir+"chi"+str(chi_num+1)+myname_alt+self.num+".xvg"
                     if os.path.exists(xvg_fn): break
               if myname == "HIS":
                  for myname_alt in ("HISA", "HISB"):
                     xvg_fn = basedir+"run"+sim_index_str+chi_dir+"chi"+str(chi_num+1)+myname_alt+self.num+".xvg"
                     if os.path.exists(xvg_fn): break

            (data,titlestr)=readxvg(xvg_fn)
            self.numangles[sequential_sim_num] = len(data[:,1])
            #print "sim_num "+str(sim_num)+" numangles: "+str(self.numangles[sim_num])
            data[:,1] = (data[:,1] + 180)%360 - 180 #Check and make sure within -180 to 180; I think this should do it 
            self.angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = data[:,1]
      ### Rank-order data for adaptive partitioning -- note that I give the rank over the pooled data for each angle in each sim
      shifted_angles = self.angles.copy()
      shifted_angles[shifted_angles < 0] += 360 # wrap these around so that -180 and 180 are part of the same bin
      for chi_num in range(self.nchi):
         #print "Chi:"+str(mychi+1)+"\n"
         shifted_flat = resize(((swapaxes((shifted_angles[chi_num,:,:]).copy(),0,1))),(sum(self.numangles)))
         sorted_angles = sort(shifted_flat)
         #print "numangles: "+str(self.numangles)+" sum_numangles: "+str(sum(self.numangles))+" num of sorted angles: "+str(len(sorted_angles))+" num of datapoints: "+str(len(shifted_flat))
         #print (searchsorted(sorted_angles,shifted_angles[chi_num,0,:]))
         for sequential_sim_num in range(num_sims):
             #print "Shifted Angles: Length: "+str(shape(shifted_angles[chi_num,0,:]))+"\n"
             #print shifted_angles[chi_num,sim_num,:]
             #print "Sorted Angles: Length\n"+str(shape(sorted_angles))+"\n"
             #print sorted_angles[:]
             #print "Rank ordered Angles \n"
             self.rank_order_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
                                      searchsorted(sorted_angles,shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles
             #print self.rank_order_angles[chi_num,sim_num,:]
	     #print "\n"
   
   # load angle info from an all_angle_info object
   def _load_pdb_data(self, all_angle_info, max_angles):
      #shifted_angles = zeros((self.nchi,all_angle_info.num_sims,max_angles),float32)
      #shifted_angles[:,:,:] = -999 #a null value other than zero      
      for sequential_sim_num in range(self.num_sims): #range(all_angle_info.num_sims):
          curr_angles, numangles = all_angle_info.get_angles(sequential_sim_num, self.sequential_num, self.name, self.backbone_only, self.phipsi)
          self.angles[:, sequential_sim_num, 0:numangles] = (curr_angles + 180)%360 - 180
          self.numangles[sequential_sim_num]= numangles
      ### Rank-order data for adaptive partitioning -- note that I give the rank over the pooled data for each angle in each sim
      shifted_angles = self.angles.copy()
      shifted_angles[shifted_angles < -90] += 360 # wrap these around so that -180 and 180 are part of the same bin
      for chi_num in range(self.nchi):
         #print "Chi:"+str(mychi+1)+"\n"
         shifted_flat = resize(((swapaxes((shifted_angles[chi_num,:,:]).copy(),0,1))),(sum(self.numangles)))
         sorted_angles = sort(shifted_flat)
         #print "numangles: "+str(self.numangles)+" sum_numangles: "+str(sum(self.numangles))+" num of sorted angles: "+str(len(sorted_angles))+" num of datapoints: "+str(len(shifted_flat))
         #print shifted_angles[chi_num,0,:200]
         #print sorted_angles[:200]
         #print (searchsorted(sorted_angles,shifted_angles[chi_num,0,:]))
         for sequential_sim_num in range(self.num_sims):
             self.rank_order_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
                                   searchsorted(sorted_angles,shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles
             #print self.rank_order_angles[chi_num,sim_num,:]
            
   # load the angles and calculate the total entropy of this residue and its variance.
   # The entropy of a residue with two chi angles is calculated as
   #     H(X1,X2) = H(X1) + H(X2) - I(X1,X2)
   # where I() is the mutual information between X1 and X2.
   # For residues with more than 3 chi angles, we assume the 2nd order approximation is sufficient.
   # The variance residue's entropy is the sum of the variances of all the terms. 
   def __init__(self,myname,mynum,sequential_num,basedir,num_sims,max_angles,xvgorpdb,binwidth,sigma=1,
                permutations=0,phipsi=0,backbone_only=0,adaptive_partitioning=0,which_runs=None,bootstrap_choose=0,
                calc_variance=True, all_angle_info=None, xvg_chidir = "/dihedrals/g_chi/"):
      self.name = myname
      self.num = mynum
      self.sequential_num = sequential_num
      self.backbone_only, self.phipsi = backbone_only, phipsi
      self.nchi=NumChis[myname] * (1 - backbone_only) + phipsi
      self.numangles = zeros((num_sims),int16)
      self.num_sims = num_sims
      self.which_runs = which_runs
      self.permutations= permutations 
      self.entropy =  zeros((self.nchi), float32)
      self.var_ent =  zeros((self.nchi), float32)
      
      print "\n#### Residue: "+self.name+" "+self.num
      binwidth = float(binwidth)
      
      bins=arange(-180,180, binwidth) #  bin edges
      nbins=len(bins) # number of bins
     
      sqrt_num_sims=sqrt(num_sims)
      self.chi_pop_hist=zeros((self.nchi,nbins),float32)
      self.chi_var_pop=zeros((self.nchi,nbins),float32)
      self.angles = zeros((self.nchi,num_sims,max_angles),float32) # the dihedral angles
      self.rank_order_angles = zeros((self.nchi,num_sims,max_angles),int16) # the dihedral angles
      ### load DATA
      if(xvgorpdb == "xvg"):
         self._load_xvg_data(basedir, num_sims, max_angles, xvg_chidir)
      if(xvgorpdb == "pdb"):
         self._load_pdb_data(all_angle_info, max_angles)

      if self.name in ("GLY", "ALA"): return

      inv_binwidth = 1.0 / binwidth
      #print "binwidth:"+str(binwidth)
      if(adaptive_partitioning == 1):
          #print "numangles:"+str(sum(self.numangles))+" binwidth = "+str(sum(self.numangles)/nbins)+"\n"
          inv_binwidth_adaptive = (nbins * 1.0 /sum(self.numangles))
          #print "inv_binwidth_adaptive: "+str(inv_binwidth_adaptive)
      
      if VERBOSE >= 2: print "Angles: ", map(int, list(self.angles[0,0,0:self.numangles[0]])), "\n\n"
      
      # use the minimum number of structures of any simulation
      # CLM: number of structures per simulation is now variable
      #self.numangles = [min(self.numangles)]*num_sims

      # find the bin for each angle and the number of counts per bin
      # need to weave this loop for speed
      self.bins = zeros((self.nchi, self.permutations + 1, num_sims, max(self.numangles)), int16) # the bin for each dihedral
      self.counts=zeros((num_sims,self.nchi,MULTIPLE_1D_BINS * nbins),float32) # normalized number of counts per bin,
      counts_marginal=zeros((num_sims,self.nchi,nbins),float32) # normalized number of counts per bin, 
      counts_adaptive=zeros((num_sims,self.nchi,nbins),float32)
      for myscramble in range(self.permutations + 1):
          for mychi in range(self.nchi):
             if (myscramble == 1): #resize
                 for sequential_sim_num in range(num_sims):
                     self.bins[mychi,:,sequential_sim_num,:] = self.bins[mychi, 0, sequential_sim_num, :] #replicate data
             if(myscramble > 0):
                  for sequential_sim_num in range(num_sims):
                      random.shuffle(self.bins[mychi, myscramble, sequential_sim_num, :])
             else:
              for sequential_sim_num in range(num_sims):
                  for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                      bin_num = binsingle(self.angles[mychi,sequential_sim_num,anglenum],MULTIPLE_1D_BINS * inv_binwidth) #no adaptive paritioning here
                      self.counts[sequential_sim_num, mychi ,bin_num] +=1 #counts are for 1-D histograms for which we use naive binning
                  self.counts[sequential_sim_num, mychi, :] /= self.numangles[sequential_sim_num] # normalize
              if(adaptive_partitioning == 0):
                  for sequential_sim_num in range(num_sims):
                      for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                          bin_num = binsingle(self.angles[mychi,sequential_sim_num,anglenum],inv_binwidth) #no adaptive paritioning here
                          counts_marginal[sequential_sim_num, mychi ,bin_num] +=1 #counts are for 1-D histograms for which we use naive binning
                          self.bins[mychi, 0, sequential_sim_num, anglenum] = bin_num #use naive binning for 2-D histograms
                      counts_marginal[sequential_sim_num,mychi, :] /= self.numangles[sequential_sim_num] # normalize
                  self.chi_pop_hist[mychi,:] = average(counts_marginal[:,mychi,:],axis=0)
              else:
                  for sequential_sim_num in range(num_sims):
                      for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                          bin_num_adaptive = binsingle_adaptive(self.rank_order_angles[mychi,sequential_sim_num,anglenum],inv_binwidth_adaptive)
                          if(bin_num_adaptive < 0):
                              print "warning!!: negative bin number wrapped to bin zero"
                              bin_num_adaptive = 0
                          if(bin_num_adaptive >= nbins):
                              print "warning!!: bin number overshot nbins, wrapping to bin nbins-1"
                              bin_num_adaptive = nbins - 1
                          self.bins[mychi, 0, sequential_sim_num, anglenum] = bin_num_adaptive  # overwrite bin value, this is used for adaptive partitioning for 2-D histograms
                          counts_adaptive[sequential_sim_num, mychi , bin_num_adaptive] +=1 #counts for 1-D histograms for 2-D mutual information calculation
                      counts_adaptive[sequential_sim_num, mychi , :] /= self.numangles[sequential_sim_num] # normalize
                  self.chi_pop_hist[mychi,:] = average(counts_adaptive[:,mychi,:],axis=0) # overwrite pop_hist value, this is used for adaptive partitioning for 2-D histograms
             

         # look out for bin values less than zero
      if len(self.bins[self.bins<0]) > 0:
          print "ERROR: bin values should not be less than zero: ",
          for i in range(num_sims):
              print arr2str( array(sorted(self.angles[mychi, i, :])))
          sys.exit(1)
      #print self.bins[0:2,:,0,:]
      calc_entropy(self.counts, num_sims, self.nchi, bootstrap_choose, calc_variance=calc_variance,entropy=self.entropy,var_ent=self.var_ent)

      #self.entropy = entropy
      #self.var_ent += var_ent
      print "Total entropy (before mutinf): %.5f (s2= %.5f )" % (sum(self.entropy), sum(self.var_ent))
      print
      # now have to correct residue entropy for correlations within the residue
      #max_S = 0.
      #for mychi1 in range(self.nchi):
      #    for mychi2 in range(mychi1+1,self.nchi):
      #       print "chi1/chi2: %d/%d" % (mychi1+1,mychi2+1)
      #       mutinf_thisdof, var_mi_thisdof, S = calc_excess_mutinf(self.chi_pop_hist[mychi1,:],self.chi_pop_hist[mychi2,:],
      #                                                              self.bins[mychi1,:,:], self.bins[mychi2,:,:], num_sims, nbins,
      #                                                              self.numangles, sigma, permutations, bootstrap_choose, calc_variance=calc_variance)
      #       max_S = max([max_S,S])
      #       print "mutinf_thisdof: "+str(mychi1)+"/"+str(mychi2)
      #       print mutinf_thisdof
      #
      #      #self.entropy -= mutinf_thisdof
      #      #self.var_ent += var_mi_thisdof
      #
      #print "Entropy: %.3f (sd= %.3f ; max(S)= %.3f ) " % (self.entropy, sqrt(self.var_ent), max_S),
      #if max_S > 0.26: print "#####",
      

# Calculate the mutual information between all pairs of residues.
# The MI between a pair of residues is the sum of the MI between all combinations of res1-chi? and res2-chi?.
# The variance of the MI is the sum of the individual variances.
# Returns mut_info_res_matrix, mut_info_uncert_matrix
def calc_pair_stats(reslist, run_params):
    rp = run_params

    #initialize the mut info matrix
    mut_info_res_matrix = zeros((len(reslist),len(reslist),6,6),float32)
    mut_info_uncert_matrix = zeros((len(reslist),len(reslist),6,6),float32)
    
    #Loop over the residue list
    for res_ind1, myres1 in zip(range(len(reslist)), reslist):
       for res_ind2, myres2 in zip(range(res_ind1, len(reslist)), reslist[res_ind1:]):
        print "\n#### Working on residues %s and %s (%s and %s):" % (myres1.num, myres2.num, myres1.name, myres2.name) 
        max_S = 0.
        if(OFF_DIAG == 1):
          for mychi1 in range(myres1.nchi):
             for mychi2 in range(myres2.nchi): 
                 print "chi1/chi2: %d/%d" % (mychi1+1,mychi2+1)
                 if(res_ind1 == res_ind2 and mychi1 == mychi2):
                     mut_info_res_matrix[res_ind1, res_ind2, mychi1, mychi2] = myres1.entropy[mychi1]
                     mut_info_uncert_matrix[res_ind1, res_ind2, mychi1, mychi2] = myres1.var_ent[mychi1]
                     max_S = 0
                 elif(res_ind1 == res_ind2 and mychi1 > mychi2):
                     mut_info_res_matrix[res_ind1, res_ind2, mychi1, mychi2] = mut_info_res_matrix[res_ind1, res_ind2, mychi2, mychi1]
                     mut_info_uncert_matrix[res_ind1, res_ind2, mychi1, mychi2] = mut_info_res_matrix[res_ind1, res_ind2, mychi2, mychi1]
                     max_S = 0
                 else:
                     mutinf_thisdof, var_mi_thisdof, S = calc_excess_mutinf(myres1.chi_pop_hist[mychi1,:],myres2.chi_pop_hist[mychi2,:],
                        myres1.bins[mychi1,:,:,:], myres2.bins[mychi2,:,:,:], rp.num_sims, rp.nbins, myres1.numangles, rp.sigma, rp.permutations, rp.bootstrap_choose, calc_variance=rp.calc_variance) 

                     mut_info_res_matrix[res_ind1 , res_ind2, mychi1, mychi2] = -0.5 * mutinf_thisdof
                     mut_info_uncert_matrix[res_ind1, res_ind2, mychi1, mychi2] = var_mi_thisdof
                     max_S = max([max_S,S])
                     mut_info_res_matrix[res_ind2, res_ind1, mychi2, mychi1] = mut_info_res_matrix[res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                     #mut_info_uncert_matrix[res_ind1, res_ind2] = mut_info_uncert_matrix[res_ind1, res_ind2]
                     mut_info_uncert_matrix[res_ind2, res_ind1, mychi2, mychi1] = mut_info_uncert_matrix[res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix

        print "mutinf=%.3f (uncert=%.3f; max(S)=%.3f" % (sum((mut_info_res_matrix[res_ind1, res_ind2, : ,:]).flatten()), sum((mut_info_uncert_matrix[res_ind1, res_ind2, :, :]).flatten()), max_S),
        if max_S > 0.26: print "#####",
        print

    return mut_info_res_matrix, mut_info_uncert_matrix

def load_resfile(run_params, load_angles=True, all_angle_info=None):
    rp = run_params
    if rp.num_structs == None: rp.num_structs = 15000

    resfile=open(rp.resfile_fn,'r')
    reslines=resfile.readlines()
    resfile.close()
    reslist = []
    #sequential_num = 0 # not used, right now...
    for resline in reslines:
       if len(resline.strip()) == 0: continue
       sequential_num, res_name, res_num = resline.split()
       #sequential_num = int(sequential_num) - 1
       if load_angles: reslist.append(ResidueChis(res_name,res_num, sequential_num, rp.xvg_basedir, rp.num_sims, rp.num_structs, rp.xvgorpdb, rp.binwidth, rp.sigma,
                                                  rp.permutations, rp.phipsi, rp.backbone_only, rp.adaptive_partitioning, rp.which_runs, calc_variance=rp.calc_variance,
                                                  all_angle_info=all_angle_info, xvg_chidir=rp.xvg_chidir))
    return reslist

# Load angle data and calculate intra-residue entropies
def load_data(run_params):
    load_resfile(run_params, load_angles=False) # make sure the resfile parses correctly (but don't load angle data yet)

    ### load trajectories from pdbs
    all_angle_info = None
    if run_params.xvgorpdb == "pdb":
       trajs = [PDBlite.PDBTrajectory(traj_fn) for traj_fn in run_params.traj_fns]
       traj_lens = array([traj.parse_len() for traj in trajs])
       run_params.num_structs = min(traj_lens)
       run_params.num_res = trajs[0].parse_len()

       all_angle_info = AllAngleInfo(run_params)
       for sequential_sim_num, true_sim_num in zip(range(run_params.num_sims), run_params.which_runs):
          all_angle_info.load_angles_from_traj(sequential_sim_num, trajs[true_sim_num-1], CACHE_TO_DISK)
       print "Shape of all angle matrix: ", all_angle_info.all_chis.shape

    ### load the residue list and angle info for those residues
    ### calculate intra-residue entropies and their variances
    reslist = load_resfile(run_params, load_angles=True, all_angle_info=all_angle_info)

    print "\n--Num residues: %d--" % (len(reslist))
    if (run_params.xvgorpdb): run_params.num_structs = max(reslist[0].numangles)

    return reslist
    
#===================================================
#READ INPUT ARGUMENTS
#===================================================
try:
   import run_profile
   run_profile.fix_args()
except: pass

usage="%prog [-t traj1:traj2] [-x xvg_basedir] resfile [simulation numbers to use]  # where resfile is in the format <1-based-index> <aa type> <res num>"
parser=OptionParser(usage)
parser.add_option("-t", "--traj_fns", default=None, type="string", help="filenames to load PDB trajectories from; colon separated (e.g. fn1:fn2)")
parser.add_option("-x", "--xvg_basedir", default=None, type="string", help="basedir to look for xvg files")
parser.add_option("-s", "--sigma", default=1.0, type="float", help="factor to use for zeroing out chi distribution pairs with high uncertainty")
parser.add_option("-w", "--binwidth", default=15.0, type="float", help="width of the bins in degrees")
parser.add_option("-n", "--num_sims", default=None, type="int", help="number of simulations")
parser.add_option("-p", "--permutations", default=20, type="int", help="number of permutations for independent mutual information, for subtraction from total Mutual Information")
parser.add_option("-d", "--xvg_chidir", default = "/dihedrals/g_chi/", type ="string", help="subdirectory under xvg_basedir/run# where chi angles are stored")
parser.add_option("-a", "--adaptive", default = "yes", type ="string", help="adaptive partitioning")
parser.add_option("-b", "--backbone", default = "chi", type = "string", help="chi: just sc  phipsi: just bb  phipsichi: bb + sc")
parser.add_option("-o", "--bootstrap_set_size", default = None, type = "int", help="perform bootstrapping within this script; value is the size of the subsets to use")
(options,args)=parser.parse_args()
if len(filter(lambda x: x==None, (options.traj_fns, options.xvg_basedir))) != 1:
    parser.error("ERROR exactly one of --traj_fns or --xvg_basedir must be specified")

# Initialize
resfile_fn = args[0]
adaptive_partitioning = (options.adaptive == "yes")
phipsi = options.backbone.find("phipsi")!=-1
backbone_only = options.backbone.find("chi")==-1
bins=arange(-180,180,options.binwidth) #Compute bin edges
nbins = len(bins)

if options.traj_fns != None:
   xvgorpdb = "pdb"
   traj_fns = options.traj_fns.split(":")
   num_sims = len(traj_fns)
else:
    xvgorpdb = "xvg"
    traj_fns = None
    num_sims = options.num_sims

if(len(args) > 1):
    which_runs = map(int, args[1:])
    num_sims = len(which_runs)
else:
    assert(num_sims != None)
    which_runs = range(1,num_sims+1)

run_params = RunParameters(resfile_fn=resfile_fn, adaptive_partitioning=adaptive_partitioning, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins,
  bootstrap_set_size=options.bootstrap_set_size, sigma=options.sigma, permutations=options.permutations, num_sims=num_sims, num_structs = None,
  binwidth=options.binwidth, bins=bins, which_runs=which_runs, xvgorpdb=xvgorpdb, traj_fns=traj_fns, xvg_basedir=options.xvg_basedir, bootstrap_choose=0, calc_variance=False, xvg_chidir=options.xvg_chidir)
print run_params

print "Calculating Entropy and Mutual Information"
  
#====================================================
#DO ANALYSIS
#===================================================

independent_mutinf_thisdof = zeros((run_params.permutations,1),float32)
timer = Timer()

### load angle data, calculate entropies and mutual informations between residues (and error-propagated variances)
if run_params.bootstrap_set_size == None:
    print run_params
    reslist = load_data(run_params)
    print "TIME to load trajectories & calculate intra-residue entropies: ", timer
    mut_info_res_matrix, mut_info_uncert_matrix = calc_pair_stats(reslist, run_params)

    prefix = run_params.get_logfile_prefix() + "_sims" + ",".join(map(str, sorted(which_runs)))
else:
    runs_superset, set_size = run_params.which_runs, run_params.bootstrap_set_size
    if set_size > len(runs_superset) or len(runs_superset) <= 1:
        print "FATAL ERROR: invalid values for bootstrap set size '%d' from runs '%s'" % (set_size, runs_superset)
        sys.exit(1)

    matrix_list = []
    run_params.calc_variance, run_params.num_sims = False, set_size
    for which_runs in xuniqueCombinations(runs_superset, set_size):
        run_params.which_runs = which_runs
        print "\n----- STARTING BOOTSTRAP RUN: %s -----" % run_params
        reslist = load_data(run_params)
        matrix_list += [calc_pair_stats(reslist, run_params)[0]] # add the entropy/mut_inf matrix to the list of matrices

    # create a master matrix
    bootstraps_mut_inf_res_matrix = zeros(list(matrix_list[0].shape) + [len(matrix_list)], float32)
    for i in range(len(matrix_list)): bootstraps_mut_inf_res_matrix[:,:,i] = matrix_list[i]
    
    mut_info_res_matrix, mut_info_uncert_matrix = bootstraps_mut_inf_res_matrix.mean(axis=2), bootstraps_mut_inf_res_matrix.std(axis=2)
    prefix = run_params.get_logfile_prefix() + "_sims%s_choose%d" % (",".join(map(str, sorted(runs_superset))), set_size)

### output results to disk
name_num_list=[]
for res in reslist: name_num_list.append(res.name + str(res.num))

if(OFF_DIAG == 1):
    output_matrix(prefix+"_mutinf.txt",mut_info_res_matrix,name_num_list,name_num_list)
    output_matrix(prefix+"_mutinf_0diag.txt",mut_info_res_matrix,name_num_list,name_num_list, zero_diag=True)
    output_matrix(prefix+"_mutinf_uncert.txt",mut_info_uncert_matrix,name_num_list,name_num_list)
if(OUTPUT_DIAG == 1):
    output_diag(prefix+"_entropy_diag_.txt",mut_info_res_matrix,name_num_list)
    output_diag(prefix+"_entropy_diag_uncert.txt",mut_info_uncert_matrix,name_num_list)

print "TIME at finish: ", timer
