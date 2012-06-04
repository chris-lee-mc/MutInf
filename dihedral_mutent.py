# MutInf program 
# Copyright 2010 Christopher McClendon and Gregory Friedland
# Released under GNU Lesser Public License

from numpy import *
import re, os, sys, os.path, time, shelve
from optparse import OptionParser
from scipy import weave
from scipy import stats as stats
from scipy import special as special
from scipy import integrate as integrate
from scipy import misc as misc
from scipy.weave import converters
import PDBlite, utils
try:
       import MDAnalysis 
except: pass

set_printoptions(linewidth=120)

mycompiler = 'gcc'  #initialize options at the global level

########################################################
##  INSTRUCTIONS:
##  
##
############################################################


##  INSTALLATION:
##  You will need to have NumPy and SciPy installed.
##  I recommend building LAPACK then ATLAS using gcc and gfortran
##  with the -fPIC option. These will provide optimized linear algebra
##  routines for your particular architecture. These plus the "weave" statements with pointer arithmetic yield performance that can approach that of C or FORTRAN code while allowing the code to be written for the most part in higher-level expressions that read more like math. 
##
##  Note: by default this code uses Numpy built with the IntelEM64-T C-compiler (command-line option "-g intelem").
##  This will use 'icc' to compile optimized inline C-code. 
##  If you instead would like to use 'gcc', indicate this with command-line option "-g gcc" 
##  PDBlite.py and utils.py should reside in the same directory as
##  dihedral_mutent.py

##  DIHEDRAL (TORSION) ANGLES FROM SIMULATIONS
##  This program uses an undersampling correction based on correlations
##  between torsions in different simulations or different equal-sized blocks of the same
##  long simulation.
##  Ideally, six or more simulations or simulation blocks are ideal.
##  This program currently assumes that all simulations or blocks are of equal length. (While in principle the code is general enough to hanlde differences in lengths of simulations/blocks, this feature is not currently enabled as a bug manifests itself when multiple sub-samples were used (i.e. when the "-o" argument is less than the "-n" argument)
##
##
##  Next, you will need dihedral angle files for each torsion, containing two space-delimited fields:
#   1) time or timestep or any number, and 2) dihedral angle in degrees.
##  These files can be named .xvg or .xvg.gz, for each torsion for each residue. for example:
##
## chi1PHE42.xvg.gz
## chi2PHE42.xvg.gz
## phiPHE42.xvg.gz
## psiPHE42.xvg.gz
##
##  I use the GROMACS utility "g_chi" to produce these files. 
##

##
## RESIDUE LIST:
## You then need a residue list file (i.e. test_new3_adaptive.reslist) with three space-delimited fields:
## 1: torsion number (i.e. 42 for above) 2. residue name  3. Number Chain (i.e. 42A for residue number 42, chain A)
## 
##
## For proteins with more than one chain, the torsion numbers can just be sequential (i.e. add 100 for each new chain), but must be unique.
##
## The residue list is used to map torsion numbers to biologically meaningful strings (i.e. "42A")
##
##  You can create a residue list using a shell command like:
## cat ${mypdb}.pdb | grep CA | awk '{print NR, substr($0,18,3), substr($0,23,4) substr($0,22,1)}' > temp

##
##  DIRECTORY STRUCTURE THE PROGRAM LOOKS FOR:
##  The -x option indicates the base directory for your dihedral data for each run.
##  Under this directory you should have directories run1, run2, run3, run4, run5, run6, for example.
##  Place the .xvg files under directories named /dihedrals/g_chi/ under each run, unless the -d option is used. For example, if "-d /" is used then the torsion angle files would be in the run1 to run6 directories.
##
##
##  RUNNING THE PROGRAM:
##  run the program like this:
##  python ~/trajtools/dihedral_mutent.py -x ~/IL2_Flex/${pdb}_62GLU_SS/ -a "yes" -o 6 -w 30 -p 0 -c "no" -n 6  test_new3_adaptive.reslist > test_new3_adaptive_out.txt
##  More options are listed in the __main__ part of the code.
##  For the purposes of this excercise, we will just aggregate all the data together as per McClendon, Friedland, et. al. 2009. 
##  To do this, the -o option (number of runs to take at a time)  needs to be the same as the -n option (the number of runs)
##  
## For this example, the directory structure would look like:
## ~/IL2_Flex/${pdb}_62GLU_SS/
##      run1/  
##           dihedrals/g_chi/*.xvg.gz
##      run2/
##           dihedrals/g_chi/*.xvg.gz
##      run3/
##           dihedrals/g_chi/*.xvg.gz
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
##
















####   CONSTANTS   #####################################################################################################################

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
OUTPUT_DIAG = 1  #set to "1" to output flat files for diagonal matrix elements
EACH_BOOTSTRAP_MATRIX = 1 #set to "1" to output flat files for each bootstrapped matrix
MULT_1D_BINS = 3 #set to "3" or more for MULT_1D_BINS times the number of 1-D bins for 1-D entropy only; single and double already done by defualt
SAFE_1D_BINS = 18 # now 20 degrees # can be approximately 2PI bins so that continuous h(x) = H(x, discrete) + log (binwidth)
FEWER_COR_BTW_BINS = 0.5 #set to "0.5" for half as many 2-D bins for 2-D correlation between tors in different sims, nbins must be an even number for this to be appropriate!
MAX_NEAREST_NEIGHBORS = 1 # up to k nearest neighbors
#SPEED = "vfast" # set to "vfast" to use the vfast_cov(); "fast" to use fast_cov(); "*" to use cov()
#print "Running at speed: ", SPEED
#################################################################################################################################
## Needs fixing if it is to be used
## Use Abel 2009 JCTC eq. 27 weighting for KNN estimates for k=1,2,3
#QKNN = zeros((3),float64)
#WEIGHTS_KNN = zeros((3),float64)
# Qk = sum_{j=k}^{inf} 1/(j^2)
#for i in range(3):
#    for j in range(999):
#        QKNN[i] += 1/((i+1+j)*(i+1+j))
#WEIGHTS_KNN = (1/QKNN) / sum(1/QKNN)


#########################################################################################################################################


#########################################################################################################################################
#### Class for two-variable vonMises distribution, used for comparing calculated/analytical entropies, mutual information ###############
#########################################################################################################################################

class vonMises:
    C = 0
    entropy = 0
    entropy_marginal1 = 0
    entropy_marginal2 = 0
    mutinf = 0
    u1 = u2 = k1 = k2 = l1 = l2 = 0
    def __k1_dC_dk1__(self):
            value = 0
            for i in range(100):
                value += misc.comb(2*i,i,exact=1) * ((self.lambd ** 2)/(4*self.k1*self.k2) ** i) * special.iv(i+1,self.k1) * special.iv(i, self.k2)
            return value * 4 * PI * PI * self.k1
    def __k2_dC_dk2__(self):
            value = 0
            for i in range(100):
                value += misc.comb(2*i,i,exact=1) * ((self.lambd ** 2)/(4*self.k1*self.k2) ** i) * special.iv(i+1,self.k2) * special.iv(i, self.k1)
            return value * 4 * PI * PI * self.k2
    def __lambd_dC_dlambd__(self):
            value = 0
            for i in range(100):
                value += misc.comb(2*i,i,exact=1) * ((self.lambd ** 2)/(4*self.k1*self.k2) ** i) * special.iv(i,self.k1) * special.iv(i, self.k2)
            return value * 4 * PI * PI
    def p(self,phi1,phi2):
        return (1.0 / self.C) * exp(self.k1 * cos(self.l1 * (phi1 - self.u1)) +
                                    self.k2 * cos(self.l2 * (phi2 - self.u2)) +
                        self.lambd * sin(self.l1 * (phi1 - self.u1)) * sin(self.l2 * (phi2 - self.u2)))
    def p1(self,phi1):
        (2*PI/self.C)*special.iv(0,sqrt(k2*k2+self.lambd*self.lambd*(sin(self.l1 * (phi1 - self.u1)) ** 2))) * exp(self.k1 * cos(self.l1 * (phi1 - self.u1)))
    def p2(self,phi2):
        (2*PI/self.C)*special.iv(0,sqrt(k1*k1+self.lambd*self.lambd*(sin(self.l2 * (phi2 - self.u2)) ** 2))) * exp(self.k2 * cos(self.l2 * (phi2 - self.u2)))
    def randarray1(self,mylen):
        myoutput = stats.uniform.rvs(size=mylen, scale=2*PI)
        for i in range(mylen):
            myoutput[i]= self.p1(myoutput[i])
        return myoutput
    def randarray2(self,mylen):
        myoutput = stats.uniform.rvs(size=mylen, scale=2*PI)
        for i in range(mylen):
            myoutput[i]= self.p2(myoutput[i])
        return myoutput
    def __init__(self,u1,u2,k1,k2,l1,l2,lambd):
        self.u1 = u1
        self.u2 = u2
        self.k1 = k1
        self.k2 = k2
        self.l1 = l1
        self.l2 = l2
        self.lambd = lambd; #lambda coupling parameter
        self.C = 0 #normalization constant 
        for i in range(100): #compute C using rapdily-converging infinite series
            self.C += misc.comb(2*i,i,exact=1) * ((self.lambd ** 2)/(4*self.k1*self.k2) ** i) * special.iv(i,self.k1) * special.iv(i, self.k2)
        #calculate entropies for joint and marginal distributions, and mutual information
        self.entropy = log(self.C) - self.k1_dC_dk1 - self.k2_dC_dk2 - self.lambd_dC_dlambd
        self.entropy_marginal1 = integrate.quad(self.p1 * log(self.p1), 0, 2*PI) - log(2*PI/self.C) - self.k1_dC_dk1
        self.entropy_marginal2 = integrate.quad(-self.p2 * log(self.p2), 0, 2*PI) - log(2*PI/self.C) - self.k2_dC_dk2
        self.mutinf = self.entropy_marginal1 + self.entropy_marginal2 - self.entropy
        #Scipy documentation on scipy functions used here:
        #scipy.misc.comb(N, k, exact=0)
        #Combinations of N things taken k at a time.
        #If exact==0, then floating point precision is used, otherwise exact long integer is computed.
        #scipy.special.iv(x1, x2[, out])
        #y=iv(v,z) returns the modified Bessel function of real order v of z. If z is of real type and negative, v must be integer valued.


#########################################################################################################################################
#####  Memory and CPU information-gathering routines  ###################################################################################
#########################################################################################################################################


"""Get some of the info from /proc on Linux
from  http://www.dcs.warwick.ac.uk/~csueda/proc.py
and http://www.velocityreviews.com/forums/t587425-python-system-information.html
'In meminfo, you're probably most interested in the fields MemTotal and
MemFree (and possibly Buffers and Cached if you want to see how much
of the used memory is in fact cache rather than user programs). You
might also want to look at SwapTotal and SwapFree.'
"""



import re
re_meminfo_parser = re.compile(r'^(?P<key>\S*):\s*(?P<value>\d*)\s*kB')


def meminfo():
    """-> dict of data from meminfo (str:int).
    Values are in kilobytes.
    """
    result = dict()
    for line in open('/proc/meminfo'):
        match = re_meminfo_parser.match(line)
        if not match:
            continue  # skip lines that don't parse
        key, value = match.groups(['key', 'value'])
        result[key] = int(value)
    #close('/proc/meminfo')
    return result



def loadavg():
    """-> 5-tuple containing the following numbers in order:
     - 1-minute load average (float)
     - 5-minute load average (float)
     - 15-minute load average (float)
     - Number of threads/processes currently executing (<= number of CPUs)
       (int)
     - Number of threads/processes that exist on the system (int)
     - The PID of the most recently-created process on the system (int)
    """
    loadavgstr = open('/proc/loadavg', 'r').readline().strip()
    data = loadavgstr.split()
    avg1, avg5, avg15 = map(float, data[:3])
    threads_and_procs_running, threads_and_procs_total = \
        map(int, data[3].split('/'))
    most_recent_pid = int(data[4])
    return avg1, avg5, avg15, threads_and_procs_running, \
             threads_and_procs_total, most_recent_pid



def cpuusage():
    """-> dict of cpuid : (usertime, nicetime, systemtime, idletime)
    cpuid "cpu" means the total for all CPUs.
    cpuid "cpuN" means the value for CPU N.
    """
    wanted_records = [line for line in open('/proc/stat') if line.startswith('cpu')]
    result = {}
    for cpuline in wanted_records:
        fields = cpuline.split()[:5]
        data = map(int, fields[1:])
        result[fields[0]] = tuple(data)
    return result



 #check for free memory at least 15%
def check_for_free_mem():
      mymemory = meminfo()
      free_memory = mymemory["MemFree"] #/ float(mymemory["SwapTotal"])
      #free_memory = mymemory["MemFree"]
      #print "Free Swap: "+str(mymemory["SwapFree"])+" Free memory: "+str(mymemory["MemFree"])+" percent free: "+str(free_memory*100)+" \% \n"
      assert(free_memory > 5*1024)


#########################################################################################################################################
### Routines for Combinatorics  #########################################################################################################
#########################################################################################################################################

# from: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
def xcombinations(items, n): # akes n distinct elements from the sequence, order matters
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcombinations(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc
def xuniqueCombinations(items, n): # takes n distinct elements from the sequence, order is irrelevant
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc
def xuniqueCombinations2(items, n): # takes n distinct elements from the sequence, order is irrelevant
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations2(items[i+1:],n-1):
                yield [items[i]]+cc

def xselections(items, n): # takes n elements (not necessarily distinct) from the sequence, order matters.
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss

#########################################################################################################################################
####  Covariance calculations: different routines for different-sized arrays so allocation happens only once ############################
#########################################################################################################################################

# calculates the covariance matrix over the rows
cov_mat = None
def fast_cov(m):
   global cov_mat

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   if cov_mat == None: cov_mat = zeros((ny, ny), float64)
   cov_mat[:,:] = 0

   code = """
      // weave1   float64
      double inv_nx, cov;

      inv_nx = 1./double(nx-1); // unbiased estimated
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
                type_converters = converters.blitz, compiler = mycompiler, runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"] )
   return cov_mat


# calculates the covariance matrix over the rows
def vfast_cov(m):
   global cov_mat

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   if cov_mat == None: cov_mat = zeros((ny, ny), float64)
   cov_mat[:,:] = 0
   
   code = """
      // weave2 float64
      double inv_nx, cov;
      double *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./double(nx-1); // unbiased estimated
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
                compiler = mycompiler,runtime_library_dirs="/usr/lib64/", library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat

#cov_mat_for1D = None
def vfast_cov_for1D(m):
   global cov_mat_for1D

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   cov_mat_for1D = zeros((ny, ny), float64)
   cov_mat_for1D[:,:] = 0
   
   code = """
      // weave3  float64
      double inv_nx, cov;
      double *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./double(nx-1); // unbiased estimated
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
                compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat_for1D

cov_mat_for1D_ind = None
def vfast_cov_for1D_ind(m):
   global cov_mat_for1D_ind

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   if cov_mat_for1D_ind == None: cov_mat_for1D_ind = zeros((ny, ny), float64)
   cov_mat_for1D_ind[:,:] = 0
   
   code = """
      // weave4 float64
      double inv_nx, cov;
      double *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./double(nx-1); // unbiased estimated
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
                compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat_for1D_ind

cov_mat_for1D_boot = None
def vfast_cov_for1D_boot(m):
   global cov_mat_for1D_boot

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   if cov_mat_for1D_boot == None: cov_mat_for1D_boot = zeros((ny, ny), float64)
   cov_mat_for1D_boot[:,:] = 0
   
   code = """
      // weave5 float64
      double inv_nx, cov;
      double *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./double(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         pCovy1 = cov_mat_for1D_boot+y1*ny;
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
            //cov_mat_for1D_boot(y1,y2) = cov*inv_nx;
            //cov_mat_for1D_boot(y2,y1) = cov_mat_for1D_boot(y1,y2);
            *(pCovy1+y2) = cov*inv_nx;
            *(cov_mat_for1D_boot+y2*ny+y1) = *(pCovy1+y2);
         }
      }
   """
   weave.inline(code, ['X', 'cov_mat_for1D_boot', 'nx', 'ny'],
                #type_converters = converters.blitz,
                compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat_for1D_boot




def vfast_cov_for1D_boot_multinomial(m):
   

   nx, ny = m.shape
   X = array(m, ndmin=2, dtype=float64)
   X -= X.mean(axis=0)

   cov_mat_for1D_boot_multinomial = zeros((ny, ny), float64)
   cov_mat_for1D_boot_multinomial[:,:] = 0
   
   code = """
      // weave5 float64
      double inv_nx, cov;
      double *pCovy1, *pXy1, *pXy2;
      int tmp;

      inv_nx = 1./double(nx-1); // unbiased estimated
      for (int y1=0; y1<ny; ++y1) {
         pCovy1 = cov_mat_for1D_boot_multinomial+y1*ny;
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
            //cov_mat_for1D_boot_multinomial(y1,y2) = cov*inv_nx;
            //cov_mat_for1D_boot_multinomial(y2,y1) = cov_mat_for1D_boot(y1,y2);
            *(pCovy1+y2) = cov*inv_nx;
            *(cov_mat_for1D_boot_multinomial+y2*ny+y1) = *(pCovy1+y2);
         }
      }
   """
   weave.inline(code, ['X', 'cov_mat_for1D_boot_multinomial', 'nx', 'ny'],
                #type_converters = converters.blitz,
                compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat_for1D_boot_multinomial


#note: here, the pairs of runs  are in the first dimension, and the permutations in the second dimension
#need to fix this routine here
cov_mat_for1Dpairs_ind = None
def vfast_var_for1Dpairs_ind(m):
   global cov_mat_for1Dpairs_ind

   nboot, nx, ny = m.shape
   #print "matrix for var: "+str(nx)+" x "+str(ny)
   X = array(m, ndmin=3, dtype=float64)
   #print "shape of matrix for var: "+str(shape(X))
   X = swapaxes(X,0,2)
   #print mean(X, axis=0)
   X -= mean(X,axis=0)
   X = swapaxes(X,0,2)

   if cov_mat_for1Dpairs_ind == None: cov_mat_for1Dpairs_ind = zeros((nboot,nx), float64)
   cov_mat_for1Dpairs_ind[:,:] = 0
   
   code = """
      // weave6 float64
      double inv_ny, var;
      double *pXy1;
      int tmp;

      inv_ny = 1./double(ny-1-1); // unbiased estimated, extra -1 to remove orig data, just want var of permuted data
      for(int bootstrap=0; bootstrap < nboot; bootstrap++) {
         for (int x1=0; x1<nx; ++x1) {
            pXy1 = X+bootstrap*ny*nx + ny*x1;
         
            var=0;
         
            for (int y=0; y<ny-1; ++y) {
               //var += X(x1,y) * X(x1,y);
               //var += *(X+x1*ny+y) * *(X+x1*ny+y);               
               //tmp=x*ny;
               var += *(pXy1+y) * *(pXy1+y);
            }
            *(cov_mat_for1Dpairs_ind + bootstrap*nx + x1) = var*inv_ny;
         
      }
     }
   """
   weave.inline(code, ['X', 'cov_mat_for1Dpairs_ind', 'nx', 'ny', 'nboot'],
                #type_converters = converters.blitz,
                compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
   return cov_mat_for1Dpairs_ind


#########################################################################################################################################
##### RunParameters class:  Stores the parameters for an Entropy/Mutinf calculation run #################################################
#########################################################################################################################################

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

#########################################################################################################################################
##### AllAngleInfo class: stores dihedral angles from pdb-format trajectories  ##########################################################
#########################################################################################################################################

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
                  print "Loaded '%s' from the cache" % key, utils.flush()
              except Exception,e:
                  print e
                  pass # if the data isn't found in the cache, then lod it

      # if it wasn't found in the cache, load it
      if not found_in_cache:
          print "Loading trajectory '%s'" % pdb_traj.traj_fn, utils.flush()

          pdb_num = 0
          for pdb in pdb_traj.get_next_pdb():
              if pdb_num >= self.num_structs: continue

              pdb_chis = pdb.calc_chis()
              pdb_phipsiomega = pdb.calc_phi_psi_omega()
              if (pdb_num+1) % 100 == 0:
                  print "Loaded pdb #%d" % (pdb_num), utils.flush()

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

   def get_angles(self, sequential_sim_num, sequential_res_num, res_name, backbone_only, phipsi, max_num_chis):
       chi_nums = []
       if phipsi == 2: chi_nums += [0,1] # which angles to use
       if backbone_only == 0: chi_nums += range(2, min(max_num_chis,NumChis[res_name]) + 2)
       #print chi_nums
       #print self.all_chis

       curr_angles = zeros((len(chi_nums), self.all_chis.shape[1]))
       for sequential_chi_num, chi_num in zip(range(len(chi_nums)), chi_nums):
           curr_angles[sequential_chi_num, :] = self.all_chis[sequential_sim_num, :, int(sequential_res_num)-1, chi_num]
       return curr_angles, curr_angles.shape[1]

#########################################################################################################################################
##### Utility Functions: Various funtions for reading input, writing output, etc.  ##########################################################
#########################################################################################################################################

#Function definition to read an xvg file
def readxvg(filename,skip,skip_over_steps):
   """Read desired xvg file; strip headers and return data as array. First column of array is times of data points; remaining columns are the data. Should properly truncate the end of the data file if any of the lines are incomplete.
INPUT: Name or path of xvg file to read.
RETURN: As a tuple:
(1) An LxN data array containing the data from the file, less the header and any aberrant lines from the end (aberrant in the sense of truncated or not following the pattern of the rest of the lines). N is the number of columns in the xvg file, and L the number of lines. It is up to the user to interpret these.
Note that units are as in xvg file (normall kJ/mol for energies from GROMACS)
(2) The title as read from the xvg file
"""

   #Read input data
   print "filename: " + filename + "\n";
   if filename.endswith(".xvg.gz"):
       import gzip
       fil = gzip.GzipFile(filename, 'r')
   elif filename.endswith(".xvg"):
       fil = open(filename,'r');
   else:
       print "ERROR: Expected and .xvg or .xvg.gz file type: " + filename
       sys.exit(1)
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
           if(len(tmp) > 3):
                  title=tmp[2]+' '+tmp[3]
        #Go to next line
        linenum+=1
   #slice off headers
   inlines=inlines[linenum:]
   #print inlines[:10] #print first 10 lines to check
   #print "\n"
   #print inlines[9990:]
   #Detect how many fields on each line in body of xvg file. 
   numfields=len(inlines[0].split())

   #Length (including any aberrant lines at the end) 
   inlength=len(inlines)
   print "inlength:"+str(inlength)+"\n"
   #Array to store data
   extra_record = 0
   if(skip == 1):
       extra_record = 0
   else:
       extra_record = 1
   dataarray=zeros((int((inlength-skip_over_steps)/skip) + extra_record,numfields),float64) #could start with zero, so add + extra_record
   
   skiplines=0
   #Read data into array
   for i in range(int((inlength-skip_over_steps)/skip) + extra_record): #could start with zero, so add + extra_record ...
      if(i*skip + skip_over_steps < inlength): # ... but make sure we don't overshoot
          entries=inlines[i*skip+skip_over_steps].split()
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
   dataarray=resize(dataarray,(int((inlength-skip_over_steps)/skip + extra_record)-skiplines,numfields))
   print "shape of dataarray:"+str(dataarray.shape)
   return (dataarray,title)


def bintouple(angle1,angle2,binwidth):
   bin1 = int(floor((angle1-0.00001 + 180) / binwidth))
   bin2 = int(floor((angle2-0.00001 + 180) / binwidth))
   return [bin1, bin2]

def binsingle(angle,inv_binwidth):
   if angle < 0: angle = 0.00000011
   if angle > 360: angle -= 360
   return int(floor((angle-0.0000001)*inv_binwidth)) #so we don't get an overshoot if angle is exactly 180

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
            #print row_num, col_num, mymatrix[row_num,col_num,row_chi,col_chi]
            myfile.write(str(mymatrix[row_num,col_num,row_chi,col_chi]))
         myfile.write(" ")
      myfile.write("\n")
   myfile.close()


def output_timeseries_chis(myfilename_prefix,myreslist,colnames, nsims = 6, nchi=6, ):
   #print "shape of matrix to ouput:"+str(mymatrix.shape)+"\n"
   min_num_angles = min(myreslist[0].numangles)
   timeseries_chis_matrix = zeros((nsims, len(reslist) * nchi, min_num_angles), float64) #initialize
   #self.angles = zeros((self.nchi,num_sims,max_angles),float64)         # the dihedral angles

   for res_ind1, myres1 in zip(range(len(reslist)), reslist):
                       print "\n#### Working on residue %s (%s):" % (myres1.num, myres1.name) , utils.flush()
                       for mychi1 in range(myres1.nchi):
                              #print "%s chi: %d/%d" % (myres1.name,int(myres1.num),mychi1+1)
                              #print "res_ind1: "+str(res_ind1)
                              #rint "mychi1: "+str(mychi1)
                              #print "nchi: " +str(nchi)
                              #print "min_num_angles: "+str(min_num_angles)
                              #print "res_ind1 * nchi + mychi1: "+str(res_ind1 * nchi + mychi1)
                              #print "myres1.angles: "
                              #print myres1.angles
                              #print "angle entries: "
                              #print myres1.angles[mychi1, :, :min_num_angles]
                              timeseries_chis_matrix[:, res_ind1 * nchi + mychi1, :] = myres1.angles[mychi1, :, :min_num_angles]
   my_file_list = []                           
   for mysim in range(nsims):
          myfile = open(myfilename_prefix + "_" + str(mysim) + ".txt",'w')
          for col_num, col_name in zip(range(len(colnames)), colnames):
                 for col_chi in range(nchi):
                        myfile.write(col_name + "_" +str(col_chi) + " ")
                 myfile.write("\n")
   
          for myrow in range(min_num_angles):
                 for col_num, col_name in zip(range(len(colnames)), colnames):
                        for col_chi in range(nchi):  
                               myfile.write(str(timeseries_chis_matrix[mysim,col_num * nchi + col_chi, myrow]))
                               myfile.write(" ")
                 myfile.write("\n")
          myfile.close()

   return timeseries_chis_matrix


def output_matrix_chis_2dhists(myfilename,mymatrix,rownames,colnames, nchi=6, nbins = 12, zero_diag=False):
    myfile = open(myfilename,'w')
    print "shape of matrix to ouput:"+str(mymatrix.shape)+"\n"
    for col_num, col_name in zip(range(len(colnames)), colnames):
        for col_chi in range(nchi):
                            for bin_j in range(nbins):
                                myfile.write(col_name + "_" +str(col_chi) +  "_" + str(bin_j) + " ")
    myfile.write("\n")
    for row_num, row_name in zip(range(len(rownames)), rownames):
        for row_chi in range(nchi):
            for col_num, col_name in zip(range(len(colnames)), colnames):
                for col_chi in range(nchi):  
                    for bin_i in range(nbins):
                        myfile.write(row_name + "_" + str(row_chi) + "_" + str(bin_i) + " ")
                        for bin_j in range(nbins):
                            if col_num == row_num and row_chi == col_chi and zero_diag:
                                myfile.write(str(0))
                            else:
                                print row_num, col_num, row_chi, col_chi, bin_i, bin_j 
                                print mymatrix[row_num,col_num,row_chi,col_chi,bin_i,bin_j]
                                myfile.write(str(mymatrix[row_num,col_num,row_chi,col_chi,bin_i,bin_j]))
                            myfile.write(" ")
            myfile.write("\n")
    myfile.close()



def make_name_num_list(reslist):
    name_num_list=[]
    for res in reslist: name_num_list.append(res.name + str(res.num))
    return name_num_list
    
def read_matrix_chis(myfilename, nchi=6, zero_diag=False):
   rownames = []
   colnames = []
   myfile = open(myfilename,'r')
   inlines = myfile.readlines()
   #print inlines
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




def compute_CA_dist_matrix(reslist,pdb_obj):  #pdb_obj is from PDBlite.py 
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


#########################################################################################################################################
##### Calc_entropy: calculates 1-dimensional entropy  ##########################################################
#########################################################################################################################################

   
# counts has shape [bootstrap_sets x nchi x nbins*MULTIPLE_1D_BINS]
def calc_entropy(counts, nchi, numangles_bootstrap, calc_variance=False, entropy = None, var_ent = None, symmetry=None):
    #assert(sum(counts[5,2,:]) >= 0.9999 and sum(counts[5,2,:]) <= 1.0001)
    nbins = counts.shape[-1] #this will be nbins*MULTIPLE_1D_BINS
    bootstrap_sets = counts.shape[0]
    
    #this normalization converts discrete space sampled to continuous domain of integral
    #using number of total bins instead of number of nonzero bins overcorrects when not all bins populated
    #consider this as the coarse-grain entropy, like the number of minima, as opposed to the shape of the minima
    #think Boltzman's law: S = k log W 
    #normalization = sum(counts > 0, axis=-1) * 1.0 / TWOPI #number of nonzero bins divided by 2pi, foreeach boot/chi
    

    normalization = nbins / (TWOPI)
    #print "log normalization"
    #print log(normalization)
    #need to expand the numangles_bootstrap over nchi and nbins for use below
    #print "counts shape: "
    #print shape(counts)
    numangles_bootstrap_chi_vector = repeat(repeat(reshape(numangles_bootstrap.copy(),(-1, 1, 1)),nbins,axis=2),nchi,axis=1)
    #print " numangles_bootstrap_chi_vector shape: "
    #print shape(numangles_bootstrap_chi_vector)
    #now replicate array of symmtery for chis, over bootstraps 
    #symmetry_bootstraps = resize(symmetry.copy(),(bootstrap_sets,nchi))
    if(VERBOSE >= 2):
        print "Counts"
        print counts
        print "entropy elements before sum and normalization:"
        print (counts * 1.0 / numangles_bootstrap_chi_vector) * (log(numangles_bootstrap_chi_vector) - special.psi(counts + SMALL) - ((-1) ** int16(counts % 2)) / (counts + 1.0))
    #print "Counts mod 2"
    #print counts % 2
    #print "Psi"
    #print special.psi(counts+SMALL)
    #print "Correction"
    #print (-((-1) ** (int16(counts % 2))) / (counts + 1.0))
    #print "entij"
    #print ((counts * 1.0 / numangles_bootstrap_chi_vector) * (log(numangles_bootstrap_chi_vector) - special.psi(counts + SMALL) - ((-1) ** int16(counts % 2)) / (counts + 1.0)))
    
    entropy[:,:] = sum((counts * 1.0 / numangles_bootstrap_chi_vector) * (log(numangles_bootstrap_chi_vector) - special.psi(counts + SMALL) - ((-1) ** int16(counts % 2)) / (counts + 1.0)),axis=2) - log(normalization)
    #the -log(normalization) is to convert from discrete space to continuous space
    #symmetry is taken into acount as a contraction of the radial space over which integration occurs
    #the +log(symmetry_bootstraps) would be to correct for torsions with symmtery,
    #but in my case I use fewer bins for the data-rich part of symmetric torsions
    #so the symmetry correction cancels out and thus does not appear here
    
         
    return entropy, var_ent


def calc_entropy_adaptive(counts_adaptive, num_sims, nchi, bootstrap_choose=0, calc_variance=False, entropy = None, var_ent = None, binwidths=None):
    #assert(entropy != None and var_ent != None)
    #if(symmetry == None): symmetry = ones((nchi),int8)
    #print counts[5,0,:]
    #print counts[5,1,:]
    #print counts[5,2,:]   
    #print sum(counts[5,0,:])
    #assert(sum(counts[5,2,:]) >= 0.9999 and sum(counts[5,2,:]) <= 1.0001)
    #nbins = counts.shape[-1] #convert historgram to normalized pdf
    #inv_normalization = 1 / normalization
    #counts *= normalization  #normalization constant (arrayed to multiply with counts elementwise) so that integral of pi ln pi is unity over [0, 2pi], symmetry already accounted for in counts
    # note: symmetry affects range and nbins of nonzero value equally
    log_counts = log(counts_adaptive)
    print "Counts_Adaptive:\n"
    print shape(counts_adaptive)
    print "Log Counts_Adaptive:\n"
    print shape(log_counts)
    print "Binwidths:\n"
    print shape(binwidths)
    #print log_counts.shape
    #print binwidths.shape
    #print shape(counts_adaptive * ( log_counts) * binwidths)
    #print shape(sum(counts_adaptive * ( log_counts) * binwidths ,axis=2))
    print  counts_adaptive * ( log_counts) * binwidths
    entropy[:,:] = -sum(counts_adaptive * ( log_counts) * binwidths ,axis=2)
    
    #print -sum(av_counts * ( log_av_counts),axis=1)
    print "Entropy:\n"
    print entropy
    
    if(calc_variance):
        #deriv_vector = - (1 + log_counts)
        #print "chi: "+str(mychi)+" bin: "+str(b)+" av_counts: " +str(av_counts)+"\n"+" entropy: "+str(-av_counts * ( log_av_counts))+" unc_av: "+str(unc_counts)+" var_ent this chi: "+str((- (1 + log_av_counts)) * (- (1 + log_av_counts)) * (unc_counts * unc_counts))
        if(num_sims > 1 and calc_variance):
            for mychi in range(nchi):
                var_ent[mychi] = vfast_cov_for1D(entropy[:,mychi])
            #bins_cov = vfast_cov_for1D(counts[:,mychi,:]) #observations are in the rows, conditions(variables) are in the columns
            #print bins_cov.shape, deriv_vector.shape
            #now we will gather all terms, diagonal and cross-terms, in the variance of joint entropy using matrix multiplication instead of explicit loop
            # var = A * Vij * T(A), where A={dy/dxj}|x=mu
            #if VERBOSE >= 2:
            #    print av_counts.shape, deriv_vector.shape, bins_cov.shape, transpose(deriv_vector).shape
            #var_ent[mychi] = inner(deriv_vector,inner(bins_cov,transpose(deriv_vector)))

         
    return entropy, var_ent



# Calculate the MI between two angle distributions, x & y.
# First calculate the vector (over all bins) of means of the joint probabilities, Pij.
# Then calculate the vectors (over all bins) of Pi*Pj and Pi+Pj.
# Now MI(x,y) = sum over all bins [ Pij * log(Pij/(Pi*Pj)) ].
# The variance of MI = d*C*d'    where d={dy/dxj}|x=mu = 1 + log(Pij/(Pi*Pj)) - (Pi+Pj) * Pij/(Pi*Pj)
#     and C is the matrix of the covariances Cij = sum over angles i [ (xi-mu(x)) * (yi-mu(y)) ]/(N-1)
pop_matrix = U = logU = pxipyj_flat = pop_matrix_sequential = U_sequential = logU_sequential = pxipyj_flat_sequential = None

    

#needed for entropy of small data sets
#### WARNING: this routine is broken and has not been tested

sumstuff_lookup = None
def sumstuff(ni_vector,numangles_bootstrap,permutations):
    nbins = ni_vector.shape[-1]
    bootstrap_sets = numangles_bootstrap.shape[0]
    mysum = zeros((bootstrap_sets,permutations+1,nbins),float64)
    maxnumangles = int(max(numangles_bootstrap))
    global sumstuff_lookup
    
    
    code_create_lookup = """
    int ni, N;
    for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
      N = *(numangles_bootstrap + mybootstrap);
      for(ni = 0; ni < N; ni ++) {
        for(int j=ni+2;j<=N+2;j++) {
            *(sumstuff_lookup + mybootstrap*maxnumangles + ni) += 1.0/(double(j)) ;
          }
      }
    }
    """
    
    code_do_lookup = """
            //weave 6b
            int ni, N;
            for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                N = *(numangles_bootstrap + mybootstrap);
                for (int permut=0; permut < permutations + 1; permut++) {
                    for(int bin=0; bin < nbins; bin++) {
                      ni = *(ni_vector  +  mybootstrap*(0 + 1)*permut*nbins  +  permut*nbins  +  bin);
                        *(mysum + mybootstrap*(0 + 1)*permut*nbins  +  permut*nbins  +  bin) =
                             *(sumstuff_lookup + mybootstrap*maxnumangles + ni);
                      
                    }
                }
            }
            """
    
    if(sumstuff_lookup == None):
        sumstuff_lookup = zeros((bootstrap_sets,maxnumangles),float64)
        weave.inline(code_create_lookup,['mysum','bootstrap_sets','permutations','ni_vector','nbins','numangles_bootstrap','maxnumangles','sumstuff_lookup'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])


    weave.inline(code_do_lookup,['mysum','bootstrap_sets','permutations','ni_vector','nbins','numangles_bootstrap','maxnumangles','sumstuff_lookup'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
    return mysum




#########################################################################################################################################
##### Mutual Information Under the Null Hypothesis via Monte Carlo Simulation  ##########################################################
#########################################################################################################################################

   
        
def calc_mutinf_multinomial_constrained( nbins, counts1, counts2, adaptive_partitioning):
    
    # allocate the matrices the first time only for speed
    
    
    
    #ninj_prior = 1.0 / float64(nbins*nbins) #Perks' Dirichlet prior
    #ni_prior = nbins * ninj_prior           #Perks' Dirichlet prior
    ninj_prior = 1.0                        #Uniform prior
    ni_prior = nbins * ninj_prior           #Uniform prior
    
    

    bootstrap_sets = 1000 # generate this many random matrices
    if(adaptive_partitioning == 0): bootstrap_sets = 100
    permutations_multinomial = 0     # just fake here

    
    numangles_bootstrap = ones((bootstrap_sets),int32) * sum(counts1[0,:])  #make a local version of this variable

    ref_counts1 = zeros((nbins),int32)
    ref_counts2 = zeros((nbins),int32)

    #create reference distribution
    code_create_ref = """
    //weave6a0
    int bin = 0;
    int mynumangles;
    mynumangles = *(numangles_bootstrap + 0); //assume only one real bootstrap set of data
    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
          if(bin >= nbins) bin=0;
          *(ref_counts1 + bin) += 1;
          *(ref_counts2 + bin) += 1;
          bin++;
    }
    """
    if(adaptive_partitioning != 0):
        weave.inline(code_create_ref,['nbins','ref_counts1','ref_counts2','numangles_bootstrap'], compiler=mycompiler)
    else:
        ref_counts1 = counts1.copy()
        ref_counts2 = counts2.copy()
    chi_counts1 = resize(ref_counts1,(bootstrap_sets,nbins))
    chi_counts2 = resize(ref_counts2,(bootstrap_sets,nbins))

    #print "chi_counts1 trial 1 :"
    #print chi_counts1[0,:]
    #print "chi_counts2 trial 1:"
    #print chi_counts2[0,:]
    chi_countdown1 = zeros((bootstrap_sets,nbins),int32)
    chi_countdown2 = zeros((bootstrap_sets,nbins),int32)

    chi_countdown1[:,:] = resize(chi_counts1[0,:].copy(),(bootstrap_sets,nbins)) #replicate up to bootstrap_sets
    chi_countdown2[:,:] = resize(chi_counts2[0,:].copy(),(bootstrap_sets,nbins)) #replicate up to bootstrap_sets

    print "counts to pick from without replacement"
    print chi_countdown1
    
    
    if(True):
       
       # 0 is a placeholder for permuations, which are not performed here; instead, analytical corrections are used
       #U = zeros((bootstrap_sets,0 + 1, nbins*nbins), float64)
       #logU = zeros((bootstrap_sets,0 + 1, nbins*nbins),float64)

       count_matrix_multi = zeros((bootstrap_sets, 0 + 1 , nbins*nbins), int32)
       
       ninj_flat_Bayes = zeros((bootstrap_sets,0 + 1,nbins*nbins),float64)
           
       
       
       numangles_bootstrap_matrix = zeros((bootstrap_sets,0+1,nbins*nbins),float64)
       for bootstrap in range(bootstrap_sets):
           for permut in range(permutations_multinomial+1):
               numangles_bootstrap_matrix[bootstrap,permut,:]=numangles_bootstrap[bootstrap]

    chi_counts1_vector = reshape(chi_counts1,(bootstrap_sets,0 + 1,nbins)) #no permuts for marginal distribution
    chi_counts2_vector = reshape(chi_counts2,(bootstrap_sets,0 + 1,nbins)) #no permuts for marginal distribution         
    numangles_bootstrap_vector = numangles_bootstrap_matrix[:,0,:nbins]    #no permuts for marginal distribution
    count_matrix_multi[:,:,:] = 0
    
    
        
    code_multi = """
    // weave6a
     #include <math.h>
     int bin1, bin2, bin1_found, bin2_found;
     for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
      int mynumangles = 0; // init
      for (int permut=0; permut < 0 + 1; permut++) {
          mynumangles = *(numangles_bootstrap + 0); //assume only one real bootstrap set of data
          for (int anglenum=0; anglenum< mynumangles; anglenum++) {
             bin1_found = 0;
             bin2_found = 0;
             while(bin1_found == 0) {       //sampling without replacement
                bin1 = int(drand48() * int(nbins));
                if( *(chi_countdown1 + mybootstrap*nbins + bin1) > 0) {
                  bin1_found = 1;
                  *(chi_countdown1 + mybootstrap*nbins + bin1) -= 1;
                }
             }
             while(bin2_found == 0) {      //sampling without replacement
                bin2 = int(drand48() * int(nbins));
                if( *(chi_countdown2 + mybootstrap*nbins + bin2) > 0) {
                  bin2_found = 1;
                  *(chi_countdown2 + mybootstrap*nbins + bin2) -= 1;
                }
             }
          //printf("bin1 %d bin2 %d\\n", bin1, bin2);
          *(count_matrix_multi  +  mybootstrap*(0 + 1)*nbins*nbins  +  permut*nbins*nbins  +  bin1*nbins + bin2) += 1;
             
          }
        }
       }
      """
    weave.inline(code_multi, ['numangles_bootstrap', 'nbins', 'count_matrix_multi','chi_counts1','chi_counts2','chi_countdown1','chi_countdown2','bootstrap_sets'],
                 #type_converters = converters.blitz,
                 compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])

    #print "done with count matrix setup for multinomial\n"
    for bootstrap in range(bootstrap_sets):
           my_flat = outer(chi_counts1[bootstrap] + ni_prior,chi_counts2[bootstrap] + ni_prior).flatten() 
           my_flat = resize(my_flat,(0 + 1,(my_flat.shape)[0]))
           ninj_flat_Bayes[bootstrap,:,:] = my_flat[:,:]

    if(numangles_bootstrap[0] > 0): #small sample stuff turned off for now cause it's broken
    #if(numangles_bootstrap[0] > 1000 and nbins >= 6):
        ent1_boots = sum((chi_counts1_vector * 1.0 / numangles_bootstrap_vector) * (log(numangles_bootstrap_vector) - special.psi(chi_counts1_vector + SMALL) - ((-1) ** (chi_counts1_vector % 2)) / (chi_counts1_vector + 1.0)),axis=2)

        ent2_boots = sum((chi_counts2_vector * 1.0 / numangles_bootstrap_vector) * (log(numangles_bootstrap_vector) - special.psi(chi_counts2_vector + SMALL) - ((-1) ** (chi_counts2_vector % 2)) / (chi_counts2_vector + 1.0)),axis=2)
        
        
        # MI = H(1)+H(2)-H(1,2)
        # where H(1,2) doesn't need absolute correction related to number of bins and symmetry factors because this cancels out in MI
        #print "Numangles Bootstrap Matrix:\n"+str(numangles_bootstrap_matrix)
        
        mutinf_thisdof = ent1_boots + ent2_boots  \
                     - (sum((count_matrix_multi * 1.0 /numangles_bootstrap_matrix)  \
                           * ( log(numangles_bootstrap_matrix) - \
                               special.psi(count_matrix_multi + SMALL)  \
                               - ((-1) ** (count_matrix_multi % 2)) / (count_matrix_multi + 1.0) \
                               ),axis=2))
        
    else:
        ent1_boots = sum((chi_counts1_vector + 1) * (1.0 / numangles_bootstrap_vector) * sumstuff(chi_counts1_vector,numangles_bootstrap_vector,permutations_multinomial),axis=2)

        ent2_boots = sum((chi_counts2_vector + 1) * (1.0 / numangles_bootstrap_vector) * sumstuff(chi_counts2_vector,numangles_bootstrap_vector,permutations_multinomial),axis=2)

        #print "shapes:"+str(ent1_boots)+" , "+str(ent2_boots)+" , "+str(sum((count_matrix_multi + 1)  * (1.0 / numangles_bootstrap_matrix) * sumstuff(count_matrix_multi,numangles_bootstrap_matrix,permutations_multinomial),axis=2))
        
        mutinf_thisdof = ent1_boots + ent2_boots \
                         -sum((count_matrix_multi + 1)  * (1.0 / numangles_bootstrap_matrix) * sumstuff(count_matrix_multi,numangles_bootstrap_matrix,permutations_multinomial),axis=2)
    
    
    
            
    
    
    
    #Now, since permutations==0 as we are doing Monte Carlo sampling instead of permutations, filter according to Bayesian estimate of distribution
    # of mutual information, M. Hutter and M. Zaffalon 2004 (or 2005).
    # Here, we will discard those MI values with p(I | data < I*) > 0.05.
    # Alternatively, we could use the permutation method or a more advanced monte carlo simulation
    # over a Dirichlet distribution to empirically determine the distribution of mutual information of the uniform
    # distribution.  The greater variances of the MI in nonuniform distributions suggest this approach
    # rather than a statistical test against the null hypothesis that the MI is the same as that of the uniform distribution.
    # The uniform distribution or sampling from a Dirichlet would be appropriate since we're using adaptive partitioning.

    #First, compute  ln(nij*n/(ni*nj) = logU, as we will need it and its powers shortly.
    #Here, use Perks' prior nij'' = 1/(nbins*nbins)
    
    if(True):  #then use Bayesian approach to approximate distribution of mutual information given data,prior
        count_matrix_multi_wprior = count_matrix_multi + ninj_prior
        numangles_bootstrap_matrix_wprior = numangles_bootstrap_matrix + ninj_prior*nbins*nbins
        numangles_bootstrap_vector_wprior = numangles_bootstrap_vector + ninj_prior*nbins*nbins
        Uij = (numangles_bootstrap_matrix_wprior) * (count_matrix_multi_wprior) / (ninj_flat_Bayes)
        logUij = log(Uij) # guaranteed to not have a zero denominator for non-zero prior (non-Haldane prior)

        Jij=zeros((bootstrap_sets, 0 + 1, nbins*nbins),float64)
        Jij = (count_matrix_multi_wprior / (numangles_bootstrap_matrix_wprior)) * logUij

        #you will see alot of "[:,0]" following. This means to take the 0th permutation, in case we're permuting data
    
        J = (sum(Jij, axis=-1))[:,0] #sum over bins ij
        K = (sum((count_matrix_multi_wprior / (numangles_bootstrap_matrix_wprior)) * logUij * logUij, axis=-1))[:,0] #sum over bins ij
        L = (sum((count_matrix_multi_wprior / (numangles_bootstrap_matrix_wprior)) * logUij * logUij * logUij, axis=-1))[:,0] #sum over bins ij
    
        #we will need to allocate Ji and Jj for row and column sums over matrix elemenst Jij:

        Ji=zeros((bootstrap_sets, 0 + 1, nbins),float64)
        Jj=zeros((bootstrap_sets, 0 + 1, nbins),float64)
        chi_counts_bayes_flat1 = zeros((bootstrap_sets,0 + 1, nbins*nbins),int32)
        chi_counts_bayes_flat2 = zeros((bootstrap_sets,0 + 1, nbins*nbins),int32)
    
    
        
        #repeat(chi_counts2[bootstrap] + ni_prior,permutations_multinomial +1,axis=0)
        for bootstrap in range(bootstrap_sets):
            chi_counts_matrix1 = reshape(resize(chi_counts1[bootstrap] + ni_prior, bootstrap_sets*(permutations_multinomial+1)*nbins),(bootstrap_sets,permutations_multinomial+1,nbins))
            chi_counts_matrix2 = reshape(resize(chi_counts2[bootstrap] + ni_prior, bootstrap_sets*(permutations_multinomial+1)*nbins),(bootstrap_sets,permutations_multinomial+1,nbins))
        
            #print "chi counts 1:" + str(chi_counts1[bootstrap])
            #print "chi counts 2:" + str(chi_counts2[bootstrap])
            
            #print "counts:\n" + str(count_matrix_multi[bootstrap])
            #now we need to reshape the marginal counts into a flat "matrix" compatible with count_matrix_multi
            # including counts from the prior
            #print "chi_counts2 shape:"+str(shape(chi_counts2))
            #print "chi_counts2[bootstrap]+ni_prior shape:"+str((chi_counts2[bootstrap] + ni_prior).shape)
            
            chi_counts_bayes_flat2[bootstrap,0,:] = repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0) #just replicate along fastest-varying axis, this works because nbins is the same for both i and j
            #print repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0)
            #handling the slower-varying index will be a little bit more tricky
            chi_counts_bayes_flat1[bootstrap,0,:] = (transpose(reshape(resize(chi_counts1[bootstrap] + ni_prior, nbins*nbins),(nbins,nbins)))).flatten() # replicate along fastest-varying axis, convert into a 2x2 matrix, take the transpose)
            #we will also need to calculate row and column sums Ji and Jj:
            Jij_2D_boot = reshape(Jij[bootstrap,0,:],(nbins,nbins))
            #Ji is the sum over j for row i, Jj is the sum over i for column j, j is fastest varying index
            #do Ji first
            Ji[bootstrap,0,:] = sum(Jij_2D_boot, axis=1)
            Jj[bootstrap,0,:] = sum(Jij_2D_boot, axis=0)


        #ok, now we calculated the desired quantities using the matrices we just set up
    
    
        numangles_bootstrap_wprior = numangles_bootstrap + ninj_prior*nbins*nbins
        
        M = (sum((1.0/(count_matrix_multi_wprior + SMALL) - 1.0/chi_counts_bayes_flat1 -1.0/chi_counts_bayes_flat2 \
                  + 1.0/numangles_bootstrap_matrix_wprior) \
                 * count_matrix_multi_wprior * logUij, axis=2))[:,0]

        Q = (1 - sum(count_matrix_multi_wprior * count_matrix_multi_wprior / ninj_flat_Bayes, axis=2))[:,0]

        if(VERBOSE >= 2):
            print "shapes"
            print "Ji:   "+str(shape(Ji))
            print "Jj:   "+str(shape(Jj))
            print "chi counts matrix 1:   "+str(chi_counts_matrix1)
            print "chi counts matrix 2:   "+str(chi_counts_matrix2)
            print "numangles bootstrap wprior:   "+str(shape(numangles_bootstrap_wprior))
            print "count_matrix_multi_wprior:" +str(shape(count_matrix_multi_wprior))
            print "chi_counts_bayes_flat1:"+str(shape(chi_counts_bayes_flat1))
            print "chi_counts_bayes_flat2:"+str(shape(chi_counts_bayes_flat2))
            
        P = (sum((numangles_bootstrap_vector_wprior) * Ji * Ji / chi_counts_matrix1, axis=2) \
             + sum(numangles_bootstrap_vector_wprior * Jj * Jj / chi_counts_matrix2, axis=2))[:,0]

        #Finally, we are ready to calculate moment approximations for p(I | Data)
        E_I_mat = ((count_matrix_multi_wprior) / (numangles_bootstrap_matrix_wprior * 1.0)) \
                  * (  special.psi(count_matrix_multi_wprior + 1.0) \
                     - special.psi(chi_counts_bayes_flat1 + 1.0) \
                     - special.psi(chi_counts_bayes_flat2 + 1.0) \
                     + special.psi(numangles_bootstrap_matrix_wprior + 1.0))
        E_I = average(sum(E_I_mat,axis = 2)[:,0]) # to get rid of permutations dimension and to average over bootstrap samples
        Var_I_runs = var(sum(E_I_mat,axis = 2)[:,0]) #average over bootstraps is needed here
        
        Var_I = abs(average( \
             ((K - J*J))/(numangles_bootstrap_wprior + 1) + (M + (nbins-1)*(nbins-1)*(0.5 - J) - Q)/((numangles_bootstrap_wprior + 1)*(numangles_bootstrap_wprior + 2)), axis=0)) # all that remains here is to average over bootstrap samples

        #now for higher moments, leading order terms

        E_I3 = average((1.0 / (numangles_bootstrap_wprior * numangles_bootstrap_wprior) ) \
               * (2.0 * (2 * J**3 -3*K*J + L) + 3.0 * (K + J*J - P)))

        E_I4 = average((3.0 / (numangles_bootstrap_wprior * numangles_bootstrap_wprior)) * ((K - J*J) ** 2))

        #convert to skewness and kurtosis (not excess kurtosis)

        
        print "Descriptive Mutinf Multinomial:"
        print average(mutinf_thisdof[:,0])
        #print "Moments for Bayesian p(I|Data) for Multinomial Dist, w/constrained marginal counts:"
        print "E_I      Multinomial: "+str(E_I)
        #print "Var_I    Multinomial: "+str(Var_I)
        #print "Stdev_I  Multinomial:"+str(sqrt(Var_I))
        #print "skewness Multinomial: "+str(E_I3/ (Var_I ** (3/2)) )
        #print "kurtosis Multinomial: "+str(E_I4/(Var_I ** (4/2)) )
        
        #print "mutinf multinomial shape:"
        #print mutinf_thisdof.shape
        print
    

    else:  #do nothing
        print 

    
    var_mi_thisdof = vfast_cov_for1D_boot_multinomial(reshape(mutinf_thisdof[:,0].copy(),(mutinf_thisdof.shape[0],1)))[0,0]
    
    
    
    return E_I, Var_I , E_I3, E_I4, Var_I_runs, mutinf_thisdof[:,0], var_mi_thisdof



#########################################################################################################################################
##### Mutual Information From Joint Histograms, with calculation of MI under Null Hypothesis ############################################
#########################################################################################################################################

   

count_matrix = None
count_matrix_star = None
E_I_multinomial = None
Var_I_multinomial = None
E_I3_multinomial = None
E_I4_multinomial = None
Var_I_runs_multinomial = None
mutinf_multinomial = None
var_mutinf_multinomial = None
mutinf_multinomial_sequential = None
var_mutinf_multinomial_sequential = None
E_I_uniform = None
numangles_bootstrap_matrix = None
numangles_bootstrap_vector = None
#ninj_flat = None
#ninj_flat_Bayes = None
#here the _star terms are for alternative ensemble weighting, such as when termsl like p* ln p are desired
#NOTE: it is assumed here that both the regular arrays and the star arrays have the same number of bootstraps, it is the job of the calling routine to ensure this
def calc_mutinf_corrected(chi_counts1, chi_counts2, bins1, bins2, chi_counts_sequential1, chi_counts_sequential2, bins1_sequential, bins2_sequential, num_sims, nbins, numangles_bootstrap, numangles, calc_variance=False,bootstrap_choose=0,permutations=0,which_runs=None,pair_runs=None, calc_mutinf_between_sims="yes", file_prefix=None, plot_2d_histograms=False, adaptive_partitioning = 0, chi_counts1_star=None, chi_counts2_star=None, bins1_star=None, bins2_star=None,chi_counts_sequential1_star=None, chi_counts_sequential2_star=None, bins1_sequential_star=None, bins2_sequential_star=None, numangles_star=None, numangles_bootstrap_star=None, bootstrap_choose_star=None ):
    global count_matrix, count_matrix_sequential, ninj_flat_Bayes, ninj_flat_Bayes_sequential # , ninj_flat
    global nbins_cor, min_angles_boot_pair_runs, numangles_bootstrap_matrix, numangles_boot_pair_runs_matrix, numangles_bootstrap_vector
    global min_angles_boot_pair_runs_matrix, min_angles_boot_pair_runs_vector
    global count_matrix_star, count_matrix_sequential_star, ninj_flat_Bayes_star, ninj_flat_Bayes_sequential_star, numangles_bootstrap_matrix_star
    global E_I_multinomial, Var_I_multinomial, E_I3_multinomial, E_I4_multinomial, Var_I_runs_multinomial 
    global mutinf_multinomial, var_mutinf_multinomial, mutinf_multinomial_sequential, var_mutinf_multinomial_sequential
    global E_I_uniform
    # allocate the matrices the first time only for speed
    if(bootstrap_choose == 0):
        bootstrap_choose = num_sims
    bootstrap_sets = len(which_runs)
    assert(all(chi_counts1 >= 0))
    assert(all(chi_counts2 >= 0))
    permutations_sequential = permutations
    #permutations_sequential = 0 #don't use permutations for mutinf between sims
    max_num_angles = int(max(numangles))
    if(numangles_star != None): 
        max_num_angles_star = int(max(numangles_star))
    
    num_pair_runs = pair_runs.shape[1]
    print pair_runs
    print "bootstrap_sets: "+str(bootstrap_sets)+"num_pair_runs: "+str(num_pair_runs)+"\n"
    nbins_cor = int(nbins * FEWER_COR_BTW_BINS)
    print "nbins_cor: "+str(nbins_cor)
    #print numangles_bootstrap
    #print bins1.shape
    
    #assert(bootstrap_sets == pair_runs.shape[0] == chi_counts1.shape[0] == chi_counts2.shape[0] == bins1.shape[0] == bins2.shape[0])
    #ninj_prior = 1.0 / float64(nbins*nbins) #Perks' Dirichlet prior
    #ni_prior = nbins * ninj_prior           #Perks' Dirichlet prior
    ninj_prior = 1.0                        #Uniform prior
    ni_prior = nbins * ninj_prior           #Uniform prior
    pvalue = zeros((bootstrap_sets),float64)

    #must be careful to zero discrete histograms that we'll put data in from weaves
    #print "chi counts 1 before multinomial:"
    #print chi_counts1

    #initialize if permutations > 0
    if(permutations > 0):
           mutinf_multinomial = mutinf_multinomial_sequential = var_mutinf_multinomial_sequential = 0
    #only do multinomial if not doing permutations
    if(E_I_multinomial == None and permutations == 0 ):
        E_I_multinomial, Var_I_multinomial, E_I3_multinomial, E_I4_multinomial, Var_I_runs_multinomial, \
                         mutinf_multinomial, var_mutinf_multinomial = \
                         calc_mutinf_multinomial_constrained(nbins,chi_counts1,chi_counts2,adaptive_partitioning)
    
    if count_matrix == None:    

       count_matrix = zeros((bootstrap_sets, permutations + 1 , nbins*nbins), int32)
       count_matrix_sequential = zeros((bootstrap_sets, num_pair_runs,permutations_sequential + 1 , nbins_cor*nbins_cor), int32)
       
       
           
       
    if(numangles_bootstrap_matrix == None):
        print "numangles bootstrap: "+str(numangles_bootstrap)
        numangles_bootstrap_matrix = zeros((bootstrap_sets,permutations+1,nbins*nbins),float64)
        for bootstrap in range(bootstrap_sets):
            for permut in range(permutations+1):
                numangles_bootstrap_matrix[bootstrap,permut,:]=numangles_bootstrap[bootstrap]
        numangles_bootstrap_vector = numangles_bootstrap_matrix[:,:,:nbins]
        #print "Numangles Bootstrap Vector:\n"+str(numangles_bootstrap_vector)
    
        min_angles_boot_pair_runs_matrix = zeros((bootstrap_sets,num_pair_runs,permutations_sequential + 1,nbins_cor*nbins_cor),int32)
        min_angles_boot_pair_runs_vector = zeros((bootstrap_sets,num_pair_runs,permutations_sequential + 1,nbins_cor),int32)

        for bootstrap in range(bootstrap_sets):
            for which_pair in range(num_pair_runs):
                print "run1 "+str(pair_runs[bootstrap,which_pair,0])+" run2 "+str(pair_runs[bootstrap,which_pair,1])
                print "numangles shape:" +str(numangles.shape)
                #my_flat = outer(chi_counts_sequential1[pair_runs[bootstrap,which_pair,0]],chi_counts_sequential2[pair_runs[bootstrap,which_pair,1]]).flatten()
                #my_flat = resize(my_flat,(permutations + 1,(my_flat.shape)[0])) #replicate original data over n permutations
                min_angles_boot_pair_runs_matrix[bootstrap,which_pair,:,:] = resize(array(min(numangles[pair_runs[bootstrap,which_pair,0]],numangles[pair_runs[bootstrap,which_pair,1]]),int32),(permutations_sequential + 1, nbins_cor*nbins_cor)) #replicate original data over permutations and nbins_cor*nbins_cor
                min_angles_boot_pair_runs_vector[bootstrap,which_pair,:,:] = resize(array(min(numangles[pair_runs[bootstrap,which_pair,0]],numangles[pair_runs[bootstrap,which_pair,1]]),int32),(nbins_cor)) #replicate original data over nbins_cor, no permutations for marginal dist.            
                #ninj_flat_Bayes_sequential[bootstrap,which_pair,:,:] = my_flat[:,:]
    
    count_matrix[:,:,:] = 0
    count_matrix_sequential[:,:,:,:] =0
    # 0 is a placeholder for permuations, which are not performed here; instead, analytical corrections are used
    chi_counts1_vector = reshape(chi_counts1.copy(),(bootstrap_sets,0 + 1,nbins)) #no permutations for marginal distributions...
    chi_counts2_vector = reshape(chi_counts2.copy(),(bootstrap_sets,0 + 1,nbins)) #no permutations for marginal distributions...
    chi_counts1_vector = repeat(chi_counts1_vector, permutations + 1, axis=1)     #but we need to repeat along permutations axis
    chi_counts2_vector = repeat(chi_counts2_vector, permutations + 1, axis=1)     #but we need to repeat along permutations axis
    chi_counts1_matrix = zeros((bootstrap_sets,0 + 1,nbins*nbins),int32)        #no permutations for marginal distributions 
    chi_counts2_matrix = zeros((bootstrap_sets,0 + 1,nbins*nbins),int32)        #no permutations for marginal distributions
    
    
    for bootstrap in range(bootstrap_sets): #copy here is critical to not change original arrays
        chi_counts2_matrix[bootstrap,0,:] = repeat(chi_counts2[bootstrap].copy(),nbins,axis=-1) #just replicate along fastest-varying axis, this works because nbins is the same for both i and j
        #print repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0)
        #handling the slower-varying index will be a little bit more tricky
        chi_counts1_matrix[bootstrap,0,:] = (transpose(reshape(resize(chi_counts1[bootstrap].copy(), nbins*nbins),(nbins,nbins)))).flatten() # replicate along fastest-varying axis, convert into a 2x2 matrix, take the transpose)
        
        
    code = """
    // weave6
    // bins dimensions: (permutations + 1) * bootstrap_sets * bootstrap_choose * max_num_angles
     //#include <math.h>
     for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
      int mynumangles = 0;
      for (int permut=0; permut < permutations + 1; permut++) {
          mynumangles = *(numangles_bootstrap + mybootstrap);
          for (int anglenum=0; anglenum< mynumangles; anglenum++) {
          if(mybootstrap == bootstrap_sets - 1) {
            //printf("bin12 %i \\n",(*(bins1  +  mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))*nbins +   (*(bins2 + permut*bootstrap_sets*bootstrap_choose*max_num_angles + mybootstrap*bootstrap_choose*max_num_angles  +  anglenum)));
            }
             *(count_matrix  +  mybootstrap*(permutations + 1)*nbins*nbins  +  permut*nbins*nbins  +  (*(bins1  +  mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))*nbins +   (*(bins2 + permut*bootstrap_sets*bootstrap_choose*max_num_angles + mybootstrap*bootstrap_choose*max_num_angles  +  anglenum))) += 1;
             
          }
        }
       }
      """
    if(VERBOSE >= 2): print "about to populate count_matrix"
    weave.inline(code, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins1', 'bins2', 'count_matrix','bootstrap_sets','permutations','max_num_angles','bootstrap_choose'],
                 #type_converters = converters.blitz,
                 compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])


    ### Redundant Sanity checks
    ninj_flat = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
    ninj_flat_Bayes = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
    for bootstrap in range(bootstrap_sets):
        my_flat = outer(chi_counts1[bootstrap] + 0.0 ,chi_counts2[bootstrap] + 0.0).flatten() # have to add 0.0 for outer() to work reliably
        assert(all(my_flat >= 0))
        my_flat = resize(my_flat,(permutations + 1,(my_flat.shape)[0]))
        ninj_flat[bootstrap,:,:] = my_flat[:,:]
        #now without the Bayes prior added into the marginal distribution
        my_flat_Bayes = outer(chi_counts1[bootstrap] + ni_prior,chi_counts2[bootstrap] + ni_prior).flatten() 
        my_flat_Bayes = resize(my_flat_Bayes,(0 + 1,(my_flat_Bayes.shape)[0]))
        ninj_flat_Bayes[bootstrap,:,:] = my_flat_Bayes[:,:]
        nbins_cor = int(nbins * FEWER_COR_BTW_BINS)
    assert(all(ninj_flat >= 0))
    Pij, PiPj = zeros((nbins+1, nbins+1), float64) - 1, zeros((nbins+1, nbins+1), float64) - 1
    permutation = 0
    Pij[1:,1:]  = (count_matrix[0,permutation,:]).reshape((nbins,nbins)) 
    PiPj[1:,1:] = (ninj_flat[0,permutation,:]).reshape((nbins,nbins)) / (numangles_bootstrap[0] * 1.0)

    if(VERBOSE >= 2):
        print "First Pass:"
        print "Sum Pij: "+str(sum(Pij[1:,1:]))+" Sum PiPj: "+str(sum(PiPj[1:,1:]))
        print "Marginal Pij, summed over j:\n"
        print sum(Pij[1:,1:],axis=1)
        print "Marginal PiPj, summed over j:\n"
        print sum(PiPj[1:,1:],axis=1)   
        print "Marginal Pij, summed over i:\n"
        print sum(Pij[1:,1:],axis=0)
        print "Marginal PiPj, summed over i:\n"
        print sum(PiPj[1:,1:],axis=0)
    ### end redundant sanity checks


    #print floor(sum(Pij[1:,1:],axis=1)) == floor(sum(PiPj[1:,1:],axis=1))
    assert(abs(sum(Pij[1:,1:]) - sum(PiPj[1:,1:])) < 0.0001)
    assert(all(abs(sum(Pij[1:,1:],axis=1) - sum(PiPj[1:,1:],axis=1)) < 0.0001))
    assert(all(abs(sum(Pij[1:,1:],axis=0) - sum(PiPj[1:,1:],axis=0)) < 0.0001))

    #######################################################################################################
    ## Now for the star terms, if they are present, otherwise set the same as the normal ones
    ########################################################################################################
    if(chi_counts1_star!=None):
         if count_matrix_star == None:    

             count_matrix_star= zeros((bootstrap_sets, permutations + 1 , nbins*nbins), int32)
             #count_matrix_sequential_star = zeros((bootstrap_sets, num_pair_runs,permutations_sequential + 1 , nbins_cor*nbins_cor), int32)
         
         if(numangles_bootstrap_matrix_star == None):
             numangles_bootstrap_matrix_star = zeros((bootstrap_sets,permutations+1,nbins*nbins),float64)
             for bootstrap in range(bootstrap_sets):
                 for permut in range(permutations+1):
                     numangles_bootstrap_matrix_star[bootstrap,permut,:]=numangles_bootstrap_star[bootstrap]
         
         count_matrix_star[:,:,:] = 0
         numangles_bootstrap_vector_star = numangles_bootstrap_matrix[:,:,:nbins]
         #count_matrix_sequential_star[:,:,:,:] =0
         # 0 is a placeholder for permuations, which are not performed here; instead, analytical corrections are used
         chi_counts1_vector_star = reshape(chi_counts1_star.copy(),(bootstrap_sets,0 + 1,nbins)) #no permutations for marginal distributions...
         chi_counts2_vector_star = reshape(chi_counts2_star.copy(),(bootstrap_sets,0 + 1,nbins)) #no permutations for marginal distributions...
         chi_counts1_vector_star = repeat(chi_counts1_vector_star, permutations + 1, axis=1)     #but we need to repeat along permutations axis
         chi_counts2_vector_star = repeat(chi_counts2_vector_star, permutations + 1, axis=1)     #but we need to repeat along permutations axis
         chi_counts1_matrix_star = zeros((bootstrap_sets,0 + 1,nbins*nbins),int32)        #no permutations for marginal distributions 
         chi_counts2_matrix_star = zeros((bootstrap_sets,0 + 1,nbins*nbins),int32)        #no permutations for marginal distributions
    
         for bootstrap in range(bootstrap_sets): #copy here is critical to not change original arrays
             chi_counts2_matrix_star[bootstrap,0,:] = repeat(chi_counts2_star[bootstrap].copy(),nbins,axis=-1) #just replicate along fastest-varying axis, this works because nbins is the same for both i and j
             #print repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0)
             #handling the slower-varying index will be a little bit more tricky
             chi_counts1_matrix_star[bootstrap,0,:] = (transpose(reshape(resize(chi_counts1_star[bootstrap].copy(), nbins*nbins),(nbins,nbins)))).flatten() # replicate along fastest-varying axis, convert into a 2x2 matrix, take the transpose)

        
         code_star = """
         // weave6_star
         // bins dimensions: (permutations + 1) * bootstrap_sets * bootstrap_choose_star * max_num_angles_star
          //#include <math.h>
          for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
           int mynumangles = 0;
           for (int permut=0; permut < permutations + 1; permut++) {
               mynumangles = *(numangles_bootstrap_star + mybootstrap);
               for (int anglenum=0; anglenum< mynumangles; anglenum++) {
               if(mybootstrap == bootstrap_sets - 1) {
               //  printf("bin12 %i \\n",(*(bins1_star  +  mybootstrap*bootstrap_choose_star*max_num_angles_star  +  anglenum))*nbins +   (*(bins2_star + permut*bootstrap_sets*bootstrap_choose_star*max_num_angles_star + mybootstrap*bootstrap_choose_star*max_num_angles_star  +  anglenum)));
                 }
                   *(count_matrix_star  +  mybootstrap*(permutations + 1)*nbins*nbins  +  permut*nbins*nbins  +  (*(bins1_star  +  mybootstrap*bootstrap_choose_star*max_num_angles_star  +  anglenum))*nbins +   (*(bins2_star + permut*bootstrap_sets*bootstrap_choose_star*max_num_angles_star + mybootstrap*bootstrap_choose_star*max_num_angles_star  +  anglenum))) += 1;
             
               }
             }
             }
         """
         if(VERBOSE >= 2): print "about to populate count_matrix_star. max_num_angles_star: "+str(max_num_angles_star)+" bootstrap_choose_star:"+str(bootstrap_choose_star)+" numangles_bootstrap_star: "+str(numangles_bootstrap_star)+" bootstrap sets: "+str(bootstrap_sets)
         weave.inline(code_star, ['num_sims', 'numangles_bootstrap_star', 'nbins', 'bins1_star', 'bins2_star', 'count_matrix_star','bootstrap_sets','permutations','max_num_angles_star','bootstrap_choose_star'],
                      #type_converters = converters.blitz,
                      compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])

         
         ninj_flat_star = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
         #ninj_flat_Bayes = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
         for bootstrap in range(bootstrap_sets):
             my_flat = outer(chi_counts1_star[bootstrap] + 0.0 ,chi_counts2_star[bootstrap] + 0.0).flatten() # have to add 0.0 for outer() to work reliably
             assert(all(my_flat >= 0))
             my_flat = resize(my_flat,(permutations + 1,(my_flat.shape)[0]))
             ninj_flat_star[bootstrap,:,:] = my_flat[:,:]
             #now without the Bayes prior added into the marginal distribution
             #my_flat_Bayes = outer(chi_counts1[bootstrap] + ni_prior,chi_counts2[bootstrap] + ni_prior).flatten() 
             #my_flat_Bayes = resize(my_flat_Bayes,(0 + 1,(my_flat_Bayes.shape)[0]))
             #ninj_flat_Bayes[bootstrap,:,:] = my_flat_Bayes[:,:]
             #nbins_cor = int(nbins * FEWER_COR_BTW_BINS)
         assert(all(ninj_flat_star >= 0))
    
         #print "chi counts 1 after bin filling:"
         #print chi_counts1
    
         # print out the 2D population histograms to disk
         # erase data matrix files after creation to save space; put files in sub-directories so we don't end up with too many files in a directory for linux to handle
         Pij_star, PiPj_star = zeros((nbins+1, nbins+1), float64) - 1, zeros((nbins+1, nbins+1), float64) - 1
         mypermutation = 0
         Pij_star[1:,1:]  = (count_matrix_star[0,mypermutation,:]).reshape((nbins,nbins))
         PiPj_star[1:,1:] = (ninj_flat_star[0,mypermutation,:]).reshape((nbins,nbins)) / (numangles_bootstrap_star[0] * 1.0)
         #debug
         if(VERBOSE >= 2):
             print "Second Pass:"
             print "Sum Pij_star: "+str(sum(Pij_star[1:,1:]))+" Sum PiPj_star: "+str(sum(PiPj_star[1:,1:]))
             print "Marginal Pij_star, summed over j:\n"
             print sum(Pij_star[1:,1:],axis=1)
             print "Marginal PiPj_star, summed over j:\n"
             print sum(PiPj_star[1:,1:],axis=1)   
             print "Marginal Pij_star, summed over i:\n"
             print sum(Pij_star[1:,1:],axis=0)
             print "Marginal PiPj_star, summed over i:\n"
             print sum(PiPj_star[1:,1:],axis=0)


    if(chi_counts1_star==None):
            count_matrix_star = count_matrix
            chi_counts1_vector_star = chi_counts1_vector
            chi_counts1_matrix_star = chi_counts1_matrix
            chi_counts2_vector_star = chi_counts2_vector
            chi_counts2_matrix_star = chi_counts2_matrix
            numangles_bootstrap_vector_star = numangles_bootstrap_vector
            numangles_bootstrap_matrix_star = numangles_bootstrap_matrix
        
    
    #######################################################################################################
    
    dKLtot_dKL1_dKL2 = zeros((bootstrap_sets), float64)
    
    #if(VERBOSE >= 2):
    #    print "nonzero bins"
    #    print nonzero_bins * 1.0
    #    print (count_matrix[bootstrap,0])[nonzero_bins]
    #    print (chi_counts1_matrix[bootstrap,0])[nonzero_bins]
    #    print (chi_counts2_matrix[bootstrap,0])[nonzero_bins]
        

    #######################################################################################################
    
    assert(all(chi_counts1 >= 0))
    assert(all(chi_counts2 >= 0))
    
    ##### to work around bug, in the meantime fill count_matrix again
    count_matrix[:,:,:] = 0
    weave.inline(code, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins1', 'bins2', 'count_matrix','bootstrap_sets','permutations','max_num_angles','bootstrap_choose'],
                 #type_converters = converters.blitz,
                 compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])

    
    ninj_flat = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
    ninj_flat_Bayes = zeros((bootstrap_sets,permutations + 1,nbins*nbins),float64)
    for bootstrap in range(bootstrap_sets):
        my_flat = outer(chi_counts1[bootstrap] + 0.0 ,chi_counts2[bootstrap] + 0.0).flatten() # have to add 0.0 for outer() to work reliably
        assert(all(my_flat >= 0))
        my_flat = resize(my_flat,(permutations + 1,(my_flat.shape)[0]))
        ninj_flat[bootstrap,:,:] = my_flat[:,:]
        #now without the Bayes prior added into the marginal distribution
        my_flat_Bayes = outer(chi_counts1[bootstrap] + ni_prior,chi_counts2[bootstrap] + ni_prior).flatten() 
        my_flat_Bayes = resize(my_flat_Bayes,(0 + 1,(my_flat_Bayes.shape)[0]))
        ninj_flat_Bayes[bootstrap,:,:] = my_flat_Bayes[:,:]
        nbins_cor = int(nbins * FEWER_COR_BTW_BINS)
    assert(all(ninj_flat >= 0))
    
    #print "chi counts 1 after bin filling:"
    #print chi_counts1
    
    # print out the 2D population histograms to disk
    # erase data matrix files after creation to save space; put files in sub-directories so we don't end up with too many files in a directory for linux to handle
    Pij, PiPj = zeros((nbins+1, nbins+1), float64) - 1, zeros((nbins+1, nbins+1), float64) - 1
    mypermutation = 0
    Pij[1:,1:]  = (count_matrix[0,mypermutation,:]).reshape((nbins,nbins))
    PiPj[1:,1:] = (ninj_flat[0,mypermutation,:]).reshape((nbins,nbins)) / (numangles_bootstrap[0] * 1.0)
    #sanity checks
    if(VERBOSE >= 2):
        print "Third Pass:"
        print "Sum Pij: "+str(sum(Pij[1:,1:]))+" Sum PiPj: "+str(sum(PiPj[1:,1:]))
        print "Marginal Pij, summed over j:\n"
        print sum(Pij[1:,1:],axis=1)
        print "Marginal PiPj, summed over j:\n"
        print sum(PiPj[1:,1:],axis=1)   
        print "Marginal Pij, summed over i:\n"
        print sum(Pij[1:,1:],axis=0)
        print "Marginal PiPj, summed over i:\n"
        print sum(PiPj[1:,1:],axis=0)
    #print floor(sum(Pij[1:,1:],axis=1)) == floor(sum(PiPj[1:,1:],axis=1))
    assert(abs(sum(Pij[1:,1:]) - sum(PiPj[1:,1:])) < 0.0001)
    assert(all(abs(sum(Pij[1:,1:],axis=1) - sum(PiPj[1:,1:],axis=1)) < 0.0001))
    assert(all(abs(sum(Pij[1:,1:],axis=0) - sum(PiPj[1:,1:],axis=0)) < 0.0001))
    
    #now sum over bootstraps for 
    Counts_ij = (average(count_matrix[:,mypermutation,:], axis=0)).reshape((nbins,nbins))

    
    if plot_2d_histograms and file_prefix != None:
        file_prefix = file_prefix.replace(":", "_").replace(" ", "_")
        print file_prefix
    
        Pij[1:, 0] = PiPj[1:,0] = bins #i.e. bin cut points 0 to 360, nbins in length
        Pij[0, 1:] = PiPj[0,1:] = bins #i.e. bin cut points 0 to 360, nbins in length
    
        #Pij[1:,1:]  = average(count_matrix[:,permutation,:], axis=0).reshape((nbins,nbins))
        #PiPj[1:,1:] = average(ninj_flat_Bayes[:,permutation,:], axis=0).reshape((nbins,nbins))
        Pij[1:,1:] /= sum(Pij[1:,1:])
        PiPj[1:,1:] /= sum(PiPj[1:,1:])
    
        res1_str = "_".join(file_prefix.split("_")[:2])
        dirname = "plots_of_Pij_PiPj_nsims%d_nstructs%d_p%d_a%s/%s" % (num_sims, numangles_bootstrap[0]/len(which_runs[0]), 0, adaptive_partitioning, res1_str)
        utils.mkdir_cd(dirname)
        
        open("%s_Pij.dat"%file_prefix, "w").write(utils.arr2str2(Pij, precision=8))
        open("%s_PiPj.dat"%file_prefix, "w").write(utils.arr2str2(PiPj, precision=8))
        #open("%s_Pij_div_PiPj.dat"%file_prefix, "w").write(utils.arr2str2(Pij_div_PiPj, precision=8))
        utils.run("R --no-restore --no-save --no-readline %s_Pij.dat < ~/bin/dihedral_2dhist_plots.R" % (file_prefix))
        utils.run("R --no-restore --no-save --no-readline %s_PiPj.dat < ~/bin/dihedral_2dhist_plots.R" % (file_prefix))
        #utils.run("rm -f %s_*.dat" % file_prefix)

        utils.cd("../..")
    


    ####
    #### ***  Important Note: don't need absolute correction related to number of bins and symmetry factors because this cancels out in MI
    #### as we're using the same number of bins in each calculation *****
    ####

    #if(numangles_bootstrap[0] > 0 or nbins >= 6):
    if(True): #small sample stuff turned off for now because it's broken
    #if(numangles_bootstrap[0] > 1000 and nbins >= 6):  
        ent1_boots = sum((chi_counts1_vector_star * 1.0 / numangles_bootstrap_vector_star) * (log(numangles_bootstrap_vector) - special.psi(chi_counts1_vector + SMALL) - ((-1) ** (chi_counts1_vector % 2)) / (chi_counts1_vector + 1.0)),axis=2) 
        
        ent2_boots = sum((chi_counts2_vector_star * 1.0 / numangles_bootstrap_vector_star) * (log(numangles_bootstrap_vector) - special.psi(chi_counts2_vector + SMALL) - ((-1) ** (chi_counts2_vector % 2)) / (chi_counts2_vector + 1.0)),axis=2)
        
        
        # MI = H(1)+H(2)-H(1,2)
        # where H(1,2) doesn't need absolute correction related to number of bins and symmetry factors because this cancels out in MI
        #print "Numangles Bootstrap Matrix:\n"+str(numangles_bootstrap_matrix)
        
        mutinf_thisdof = ent1_boots + ent2_boots  \
                     - sum((count_matrix_star * 1.0 /numangles_bootstrap_matrix_star)  \
                           * ( log(numangles_bootstrap_matrix) - \
                               special.psi(count_matrix + SMALL)  \
                               - ((-1) ** (count_matrix % 2)) / (count_matrix + 1.0) \
                               ),axis=2)
    else:
        ent1_boots = sum((chi_counts1_vector + 1) * (1.0 / numangles_bootstrap_vector) * sumstuff(chi_counts1_vector,numangles_bootstrap_vector,permutations),axis=2)

        ent2_boots = sum((chi_counts2_vector + 1) * (1.0 / numangles_bootstrap_vector) * sumstuff(chi_counts2_vector,numangles_bootstrap_vector,permutations),axis=2)

        if(VERBOSE >= 2):
            print "shapes:"+str(ent1_boots)+" , "+str(ent2_boots)+" , "+str(sum((count_matrix + 1)  * (1.0 / numangles_bootstrap_matrix) * sumstuff(count_matrix,numangles_bootstrap_matrix,permutations_multinomial),axis=2))
        
        mutinf_thisdof = ent1_boots + ent2_boots \
                         -sum((count_matrix + 1)  * (1.0 / numangles_bootstrap_matrix) * sumstuff(count_matrix,numangles_bootstrap_matrix,permutations),axis=2)

    
    print "Avg Descriptive Mutinf:    "+str(average(mutinf_thisdof[:,0]))
    if(permutations == 0):
        print "Avg Descriptive MI ind:    "+str(average(mutinf_multinomial,axis=0))
    else:
        print "Avg Descriptive MI ind:    "+str(average(mutinf_thisdof[:,1:]))
        print "Number of permutations:    "+str(mutinf_thisdof.shape[1] -1)
    
    #Now, if permutations==0, filter according to Bayesian estimate of distribution
    # of mutual information, M. Hutter and M. Zaffalon 2004 (or 2005).
    # Here, we will discard those MI values with p(I | data < I*) > 0.05.
    # Alternatively, we could use the permutation method or a more advanced monte carlo simulation
    # over a Dirichlet distribution to empirically determine the distribution of mutual information of the uniform
    # distribution.  The greater variances of the MI in nonuniform distributions suggest this approach
    # rather than a statistical test against the null hypothesis that the MI is the same as that of the uniform distribution.
    # The uniform distribution or sampling from a Dirichlet would be appropriate since we're using adaptive partitioning.

    #First, compute  ln(nij*n/(ni*nj) = logU, as we will need it and its powers shortly.
    #Here, use Perks' prior nij'' = 1/(nbins*nbins)
    
    
    if(permutations == 0):  #then use Bayesian approach to approximate distribution of mutual information given data,prior
        count_matrix_wprior = count_matrix + ninj_prior
        count_matrix_wprior_star = count_matrix_star + ninj_prior #for alternative ensemble, weighting, as in cross terms like p* ln p.
        numangles_bootstrap_matrix_wprior = numangles_bootstrap_matrix + ninj_prior*nbins*nbins
        numangles_bootstrap_vector_wprior = numangles_bootstrap_vector + ninj_prior*nbins*nbins
        numangles_bootstrap_matrix_wprior_star = numangles_bootstrap_matrix_star + ninj_prior*nbins*nbins
        numangles_bootstrap_vector_wprior_star = numangles_bootstrap_vector_star + ninj_prior*nbins*nbins
        Uij = (numangles_bootstrap_matrix_wprior) * (count_matrix_wprior) / (ninj_flat_Bayes)
        logUij = log(Uij) # guaranteed to not have a zero denominator for non-zero prior (non-Haldane prior)

        Jij=zeros((bootstrap_sets, permutations + 1, nbins*nbins),float64)
        Jij = (count_matrix_wprior_star / (numangles_bootstrap_matrix_wprior)) * logUij

        #you will see alot of "[:,0]" following. This means to take the 0th permutation, in case we're permuting data
    
        J = (sum(Jij, axis=-1))[:,0] #sum over bins ij
        K = (sum((count_matrix_wprior_star / (numangles_bootstrap_matrix_wprior_star)) * logUij * logUij, axis=-1))[:,0] #sum over bins ij
        L = (sum((count_matrix_wprior_star / (numangles_bootstrap_matrix_wprior_star)) * logUij * logUij * logUij, axis=-1))[:,0] #sum over bins ij
    
        #we will need to allocate Ji and Jj for row and column sums over matrix elemenst Jij:

        Ji=zeros((bootstrap_sets, permutations + 1, nbins),float64)
        Jj=zeros((bootstrap_sets, permutations + 1, nbins),float64)
        chi_counts_bayes_flat1 = zeros((bootstrap_sets,permutations + 1, nbins*nbins),int32)
        chi_counts_bayes_flat2 = zeros((bootstrap_sets,permutations + 1, nbins*nbins),int32)
    
    
        
        #repeat(chi_counts2[bootstrap] + ni_prior,permutations +1,axis=0)
        for bootstrap in range(bootstrap_sets):
            chi_counts_matrix1 = reshape(resize(chi_counts1[bootstrap] + ni_prior, bootstrap_sets*(permutations+1)*nbins),(bootstrap_sets,permutations+1,nbins))
            chi_counts_matrix2 = reshape(resize(chi_counts2[bootstrap] + ni_prior, bootstrap_sets*(permutations+1)*nbins),(bootstrap_sets,permutations+1,nbins))
        
            if(VERBOSE >= 2):
                print "chi counts 1:" + str(chi_counts1[bootstrap])
            if(VERBOSE >= 2):
                print "chi counts 2:" + str(chi_counts2[bootstrap])

            mycounts_mat = reshape(count_matrix[bootstrap],(nbins,nbins))
            if(VERBOSE >= 2):
                print "counts:\n" + str(mycounts_mat)
            assert(all(chi_counts1[bootstrap] == sum(mycounts_mat,axis=1)))
            assert(all(chi_counts2[bootstrap] == sum(mycounts_mat,axis=0)))

            #now we need to reshape the marginal counts into a flat "matrix" compatible with count_matrix
            # including counts from the prior
            if(VERBOSE >= 2):
                print "chi_counts2 shape:"+str(shape(chi_counts2))
            #print "chi_counts2[bootstrap]+ni_prior shape:"+str((chi_counts2[bootstrap] + ni_prior).shape)
            
            chi_counts_bayes_flat2[bootstrap,0,:] = repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0) #just replicate along fastest-varying axis, this works because nbins is the same for both i and j
            #print repeat(chi_counts2[bootstrap] + ni_prior,nbins,axis=0)
            #handling the slower-varying index will be a little bit more tricky
            chi_counts_bayes_flat1[bootstrap,0,:] = (transpose(reshape(resize(chi_counts1[bootstrap] + ni_prior, nbins*nbins),(nbins,nbins)))).flatten() # replicate along fastest-varying axis, convert into a 2x2 matrix, take the transpose)
            #we will also need to calculate row and column sums Ji and Jj:
            Jij_2D_boot = reshape(Jij[bootstrap,0,:],(nbins,nbins))
            #Ji is the sum over j for row i, Jj is the sum over i for column j, j is fastest varying index
            #do Ji first
            Ji[bootstrap,0,:] = sum(Jij_2D_boot, axis=1)
            Jj[bootstrap,0,:] = sum(Jij_2D_boot, axis=0)


        #ok, now we calculated the desired quantities using the matrices we just set up
    
    
        numangles_bootstrap_wprior = numangles_bootstrap + ninj_prior*nbins*nbins
        
        M = (sum((1.0/(count_matrix_wprior + SMALL) - 1.0/chi_counts_bayes_flat1 -1.0/chi_counts_bayes_flat2 \
                  + 1.0/numangles_bootstrap_matrix_wprior) \
                 * count_matrix_wprior * logUij, axis=2))[:,0]

        Q = (1 - sum(count_matrix_wprior * count_matrix_wprior / ninj_flat_Bayes, axis=2))[:,0]

        #####DEBUGGING STATEMENTS
        #print "shapes"
        #print "Ji:   "+str(shape(Ji))
        #print "Jj:   "+str(shape(Jj))
        #print "chi counts matrix 1:   "+str(chi_counts_matrix1)
        #print "chi counts matrix 2:   "+str(chi_counts_matrix2)
        #print "numangles bootstrap wprior:   "+str(shape(numangles_bootstrap_wprior))

        P = (sum((numangles_bootstrap_vector_wprior) * Ji * Ji / chi_counts_matrix1, axis=2) \
             + sum(numangles_bootstrap_vector_wprior * Jj * Jj / chi_counts_matrix2, axis=2))[:,0]

        #####DEBUGGING STATEMENTS
        #print "ni prior:\n"+str(ni_prior)
        #print "numangles bootstrap wprior:\n"+str(numangles_bootstrap_wprior)
        #print "chi_counts_bayes_flat1:\n"+str(chi_counts_bayes_flat1)
        #print "chi_counts_bayes_flat2:\n"+str(chi_counts_bayes_flat2)

        #print intermediate values

        #print "J:"+str(J)
        #print "K:"+str(K)
        #print "L:"+str(L)
        #print "M:"+str(M)
        #print "P:"+str(P)
        #print "Q:"+str(Q)

        #Finally, we are ready to calculate moment approximations for p(I | Data)
        E_I_mat = ((count_matrix_wprior_star) / (numangles_bootstrap_matrix_wprior_star * 1.0)) \
                  * (  special.psi(count_matrix_wprior + 1.0) \
                     - special.psi(chi_counts_bayes_flat1 + 1.0) \
                     - special.psi(chi_counts_bayes_flat2 + 1.0) \
                     + special.psi(numangles_bootstrap_matrix_wprior + 1.0))
        #print "\n"
        if(VERBOSE >= 2):
            print "E_I_mat:\n"
            print E_I_mat
        E_I = (sum(E_I_mat,axis = 2))[:,0] # to get rid of permutations dimension 
        #print "Stdev of E_I over bootstraps:\n"
        #print stats.std(average(sum(E_I_mat,axis = 2), axis = 1))
        #print "estimated pop std of E_I over bootstraps, dividing by sqrt n"
        #print stats.std(average(sum(E_I_mat,axis = 2), axis = 1)) / sqrt(num_sims)
        #print "Stdev_mat:\n"
        #print sqrt(((K - J*J))/(numangles_bootstrap_wprior + 1) + (M + (nbins-1)*(nbins-1)*(0.5 - J) - Q)/((numangles_bootstrap_wprior + 1)*(numangles_bootstrap_wprior + 2)))

        Var_I = abs( \
             ((K - J*J))/(numangles_bootstrap_wprior + 1) + (M + (nbins-1)*(nbins-1)*(0.5 - J) - Q)/((numangles_bootstrap_wprior + 1)*(numangles_bootstrap_wprior + 2))) # a different variance for each bootstrap sample

        #now for higher moments, leading order terms

        E_I3 = (1.0 / (numangles_bootstrap_wprior * numangles_bootstrap_wprior) ) \
               * (2.0 * (2 * J**3 -3*K*J + L) + 3.0 * (K + J*J - P))

        E_I4 = (3.0 / (numangles_bootstrap_wprior * numangles_bootstrap_wprior)) * ((K - J*J) ** 2)

        #convert to skewness and kurtosis (not excess kurtosis)

        
        
        print "Moments for Bayesian p(I|Data):"
        print "E_I:                   "+str(E_I)
        print "Var_I:                 "+str(Var_I)
        #print "Stdev_I:               "+str(sqrt(Var_I))
        #print "skewness:              "+str(E_I3/ (Var_I ** (3/2)) )
        #print "kurtosis:              "+str(E_I4/(Var_I ** (4/2)) )

        
        
        def Edgeworth_pdf(u1,u2,u3,u4):
            #convert moments to cumulants for u0=1, u1=0, u2=1, central moments normalized to zero mean and unit variance
            skewness = u3 / (u2 ** (3/2))
            excess_kurtosis = u4 / (u2 ** (4/2)) - 3
            s = sqrt(u2)
            k3 = skewness
            k4 = excess_kurtosis
            
            return lambda x: stats.norm.pdf((x-u1)/s) * (1.0 + (1.0/6)*k3*(special.hermitenorm(3)((x-u1)/s)) \
                              + (1.0/24)*k4*(special.hermitenorm(4)((x-u1)/s)))

        def Edgeworth_quantile(crit_value,u1,u2,u3,u4):
            #convert critical value to z-score
            #func = Edgeworth_pdf(u1,u2,u3,u4)
            #print func
            #normalization = integrate.quad( func, -100*sqrt(u2), 100*sqrt(u2) )[0]
            #integral = integrate.quad( func, -integrate.inf, crit_value)[0] #output just the definite integral
            #print "integral:              "+str(integral)
            #print "normalization:         "+str(normalization)
            #pval = abs(integral)
            #print "plusminus_sigma_check:"+str(integrate.quad( func, u1-sqrt(u2), u1+sqrt(u2))[0]/normalization)
            pval_gauss = stats.norm.cdf((crit_value - u1)/sqrt(u2))
            #print "p-value Bayes Edgeworth: "+str(pval)
            print "p-value Bayes Gaussian:  "+str(pval_gauss)
            return pval_gauss

        #if(bootstrap_sets > 1):
        #    numangles_bootstrap_avg_wprior = average(numangles_bootstrap_wprior)
        #else:
        #    numangles_bootstrap_avg_wprior = numangles_bootstrap_wprior[0]

        if E_I_uniform == None:
            E_I_uniform = average(special.psi(numangles_bootstrap_wprior / (nbins * nbins) + 1) \
                                  - special.psi(numangles_bootstrap_wprior / nbins + 1) \
                                  - special.psi(numangles_bootstrap_wprior / nbins + 1) \
                                  + special.psi(numangles_bootstrap_wprior + 1))

        #####DEBUGGING STATEMENTS
        #print "Edgeworth pdf for E_I and E_I_uniform:"


        
        #print Edgeworth_pdf( E_I, Var_I, E_I3, E_I4)(E_I)
        #print Edgeworth_pdf( E_I, Var_I, E_I3, E_I4)(E_I_uniform)

        
        
        #print "E_I_uniform            :"+str(E_I_uniform)
        #print "E_I_multinomial_constr :"+str(E_I_multinomial)
        #now, determine the probability given the data that the true mutual information is
        #greater than that expected for the uniform distribution plus three sigma

        #pvalue for false positive
        #lower pvalue is better
        ##############################
        
        
        for bootstrap in range(bootstrap_sets):
            #for speed, use shortcuts for obviously significant or obviously insignificant mutual information
            if (E_I[bootstrap] < E_I_multinomial):
                pvalue[bootstrap] = 1.0
            else:
                if(E_I[bootstrap] > E_I_multinomial + 10 * sqrt(Var_I[bootstrap])):
                    pvalue[bootstrap] = 0.0
                else:
                    pvalue[bootstrap] = Edgeworth_quantile(E_I_multinomial, E_I[bootstrap], Var_I[bootstrap] , E_I3[bootstrap], E_I4[bootstrap])
            print "Bayesian P(I<E[I]mult) bootstrap sample:"+str(bootstrap)+" = "+str(pvalue[bootstrap])
            num_greater_than_obs_MI = sum(1.0 * (mutinf_multinomial > mutinf_thisdof[bootstrap,0]))
            if num_greater_than_obs_MI < 1:
                num_greater_than_obs_MI = 0
            pvalue_multinomial = num_greater_than_obs_MI * 1.0 / float32(mutinf_multinomial.shape[0])
            #print "Mutinf Multinomial Shape:"+str(mutinf_multinomial.shape)
            print "Num Ind Greater than Obs:"+str(num_greater_than_obs_MI)
            print "Descriptive P(I=I_mult):"+str(pvalue_multinomial)
            pvalue[bootstrap] = max(pvalue[bootstrap], pvalue_multinomial)
            print "Max pvalue             :"+str(pvalue[bootstrap])
        if(VERBOSE >= 2):
            print "integrate check: "+str(Edgeworth_quantile(integrate.inf,  E_I, Var_I, E_I3, E_I4))
        #could use a monte carlo simulation to generate distribution of MI of uniform distribution over adaptive partitioning
        #this would be MI of independent variables
        #for non-Bayesian significance test against null hypothesis of independence
        #but not at the present time
        

    else:  #use permutation test to filter out true negatives, possibly in addition to the Bayesian filter above
        
        #pvalue is the fraction of mutinf values from samples of permuted data that are greater than the observed MI
        #pvalue for false negative
        #lower pvalue is better
        if(permutations > 0): #otherwise, keep pvalue as 0 for now, use wilcoxon signed ranks test at the end
            for bootstrap in range(bootstrap_sets):
                num_greater_than_obs_MI = sum(mutinf_thisdof[bootstrap,1:] > mutinf_thisdof[bootstrap,0])
                print "number of permutations with MI > MI(observed): "+str(num_greater_than_obs_MI)
                pvalue[bootstrap] += num_greater_than_obs_MI * 1.0 / permutations
                print "Descriptive P(avg(I) = avg(I,independent)"+str(pvalue[bootstrap])
                Var_I = 0 #will be overwritten later
        else:
            for bootstrap in range(bootstrap_sets):
                num_greater_than_obs_MI = sum(1.0 * (mutinf_multinomial > mutinf_thisdof[bootstrap,0]))
                print "number of permutations with MI > MI(observed): "+str(num_greater_than_obs_MI)
                pvalue_multinomial = num_greater_than_obs_MI * 1.0 / float32(mutinf_multinomial.shape[0])
                print "Descriptive P(I=I_mult):"+str(pvalue_multinomial)
                pvalue[bootstrap] = max(pvalue[bootstrap], pvalue_multinomial)
                print "Descriptive P(avg(I) = avg(I,independent)"+str(pvalue[bootstrap])
                Var_I = 0 #will be overwritten later

    
    #print mutinf_thisdof
    if(bootstrap_sets > 1 and calc_variance == True):
       if(permutations > 0):
           var_mi_thisdof = (vfast_cov_for1D_boot(reshape(mutinf_thisdof[:,0],(mutinf_thisdof.shape[0],1))) - average(mutinf_thisdof[:,1:],axis=1))[0,0]
       else:
           var_mi_thisdof = (vfast_cov_for1D_boot(reshape(mutinf_thisdof[:,0],(mutinf_thisdof.shape[0],1))))[0,0]
    else:
       var_mi_thisdof = sum(Var_I)
    print "var_mi_thisdof: "+str(var_mi_thisdof)+"\n"
    if(calc_mutinf_between_sims == "no" or num_pair_runs <= 1):
        mutinf_thisdof_different_sims= zeros((bootstrap_sets,permutations_sequential + 1),float64)
        return mutinf_thisdof, var_mi_thisdof , mutinf_thisdof_different_sims, 0, average(mutinf_multinomial), 0, pvalue, dKLtot_dKL1_dKL2, Counts_ij

    #########################################################################################################
    ##
    ##  Now we will calculate mutinf between torsions in between sims for the undersampling correction 
    ##
    ##
    ##  For now, we will use the same undersampling correction for both regular and "_star" terms
    ##  
    ## 
    #########################################################################################################

    count_matrix_sequential[:,:,:,:] = 0
    
    max_num_angles = int(max(numangles))

    if(permutations == 0):
        if(mutinf_multinomial_sequential == None or adaptive_partitioning == 0):
            E_I_multinomial_sequential, Var_I_multinomial_sequential, E_I3_multinomial_sequential, E_I4_multinomial_sequential, Var_I_runs_multinomial_sequential, mutinf_multinomial_sequential, var_mutinf_multinomial_sequential = \
                                        calc_mutinf_multinomial_constrained(nbins_cor,chi_counts_sequential1.copy(),chi_counts_sequential2.copy(),adaptive_partitioning)
    else:
        mutinf_multinomial_sequential = var_mutinf_multinomial_sequential = 0
    
    
    code = """
     // weave7
     #include <math.h>
     int mynumangles, mynumangles1, mynumangles2, run1, run2, fetch1, fetch2, bin1, bin2 = 0;
     for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
       for(int which_run_pair=0; which_run_pair < num_pair_runs; which_run_pair++) {
         mynumangles1 = *(numangles + (*(pair_runs + mybootstrap*num_pair_runs + which_run_pair*2 + 0)));
         mynumangles2 = *(numangles + (*(pair_runs + mybootstrap*num_pair_runs + which_run_pair*2 + 1)));
         if(mynumangles1 <= mynumangles2) mynumangles = mynumangles1; else mynumangles = mynumangles2;
         run1 = (*(pair_runs + mybootstrap*num_pair_runs + which_run_pair*2  + 0));
         run2 = (*(pair_runs + mybootstrap*num_pair_runs + which_run_pair*2 + 1));
         for (int permut=0; permut < permutations_sequential + 1; permut++) {
           for (int anglenum=0; anglenum< mynumangles; anglenum++) {
             bin1 = *(bins1_sequential  +  run1*max_num_angles  +  anglenum);
             bin2 = *(bins2_sequential  + permut*num_sims*max_num_angles + run2*max_num_angles  +  anglenum);
             *(count_matrix_sequential  +  mybootstrap*num_pair_runs*(permutations_sequential + 1)*nbins_cor*nbins_cor + which_run_pair*(permutations_sequential + 1)*nbins_cor*nbins_cor  +  permut*nbins_cor*nbins_cor  +  bin1*nbins_cor +   bin2) += 1;
          }
        }
       }
       }
    """                                           
                                          
                                           
    weave.inline(code, ['num_sims', 'numangles','max_num_angles', 'nbins_cor', 'bins1_sequential', 'bins2_sequential', 'count_matrix_sequential','pair_runs','num_pair_runs','bootstrap_sets','permutations_sequential'],
                 #type_converters = converters.blitz,
                 compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
    #for bootstrap in range(bootstrap_sets):
    #    for which_run_pair
    #    count_matrix_sequential[bootstrap,:,:] /= min((numangles[pair_runs[bootstrap,0]], numangles[pair_runs[bootstrap,1]]))
    

                                                                              
                                     
     
    
    #logU_sequential[:,:,:] = 0
    #U_sequential = count_matrix_sequential / (ninj_flat_Bayes_sequential + SMALL)
    if(VERBOSE >= 2):
        print "count_matrix_sequential\n"
        print count_matrix_sequential[0,0,0:]

    ##### DEBUGGING STUFF ##############################
    #print "chi counts sequential1"
    #print chi_counts_sequential1[0]
    #print "shape of count_matrix_sequential\n"
    #print shape(count_matrix_sequential)
    #print "chi pop hist sequential1\n"
    #print chi_pop_hist_sequential1[pair_runs[0,0,0]]*numangles[pair_runs[0,0,0]]
    #print "chi pop hist sequential1\n"
    #print chi_pop_hist_sequential1[0]*numangles[0]
    #print "sum"
    #print sum(chi_pop_hist_sequential1[pair_runs[0,0,0]]*numangles[pair_runs[0,0,0]])
    #print "chi pop hist sequential2\n"
    #print chi_pop_hist_sequential2[pair_runs[0,0,1]]*numangles[pair_runs[0,0,1]]
    #print "chi pop hist sequential2\n"
    #print chi_pop_hist_sequential2[1]*numangles[1]
    #print "sum"
    #print sum(chi_pop_hist_sequential2[1]*numangles[1])
    #print "before log\n"
    #logU_sequential = log(U_sequential + SMALL)
                                                                              
    #print "min angles boot pair runs vector shape:"
    #print  shape(min_angles_boot_pair_runs_vector)
    ######################################################


    #take entropy for sims in chi counts sequential, summing over sims in each bootstrap
    ent1_boots_sims = zeros((bootstrap_sets, num_pair_runs, permutations_sequential + 1),float64)
    ent2_boots_sims = zeros((bootstrap_sets, num_pair_runs, permutations_sequential + 1),float64)

    for bootstrap in range(bootstrap_sets):
        for which_pair in range(num_pair_runs):
            #pick index 0 for axis=2 at the end because the arrays are over 
            
            myent1                                \
                                                  =\
                                                  sum((chi_counts_sequential1[pair_runs[bootstrap,which_pair,0]] \
                                                       * 1.0 / min_angles_boot_pair_runs_vector[bootstrap,pair_runs[bootstrap,which_pair,0],0,:]) \
                                                      * (log(min_angles_boot_pair_runs_vector[bootstrap,pair_runs[bootstrap,which_pair,0],0,:]) \
                                                         - special.psi( \
                chi_counts_sequential1[pair_runs[bootstrap,which_pair,0]] + SMALL) \
                                                         - ((-1) ** ( \
                int32(chi_counts_sequential1[pair_runs[bootstrap,which_pair,0]]) % 2)) / (chi_counts_sequential1[pair_runs[bootstrap,which_pair,0]] + 1.0)),axis=0) #sum over bins
            
            #print "ent1 boots sims thispair shape: "
            #print myent1.shape
            #print "ent1 boots: "+str(myent1)

            ent1_boots_sims[bootstrap,which_pair,:] = myent1
            
            myent2                                \
                                                  =\
                                                  sum((chi_counts_sequential2[pair_runs[bootstrap,which_pair,1]] \
                                                       * 1.0 / min_angles_boot_pair_runs_vector[bootstrap,pair_runs[bootstrap,which_pair,1],0,:]) \
                                                      * (log(min_angles_boot_pair_runs_vector[bootstrap,pair_runs[bootstrap,which_pair,1],0,:]) \
                                                         - special.psi( \
                chi_counts_sequential2[pair_runs[bootstrap,which_pair,1]] + SMALL) \
                                                         - ((-1) ** ( \
                int32(chi_counts_sequential2[pair_runs[bootstrap,which_pair,1]]) % 2)) / (chi_counts_sequential2[pair_runs[bootstrap,which_pair,1]] + 1.0)),axis=0) #sum over bins


            #print "ent1 boots sims thispair shape: "
            #print myent2.shape
            #print "ent1 boots: "+str(myent2)

            ent2_boots_sims[bootstrap,which_pair,:] = myent2
            #ent1_boots = average(ent1_boots_sims,axis=1) # avg over sims in each bootstrap
            #ent2_boots = average(ent2_boots_sims,axis=1) # avg over sims in each bootstrap

    #print "ent1 boots sims shape:"
    #print shape(ent1_boots_sims)
    #print "ent1 boots sims :"
    #print ent1_boots_sims
    #print "ent2 boots sims :"
    #print ent2_boots_sims
    ent1_boots = ent1_boots_sims
    ent2_boots = ent2_boots_sims
    
    
    
    # MI = H(1)+H(2)-H(1,2)
    # where H(1,2) doesn't need absolute correction related to number of bins and symmetry factors because this cancels out in MI

    ent12_boots = sum((count_matrix_sequential * 1.0 /min_angles_boot_pair_runs_matrix)  \
                           * ( log(min_angles_boot_pair_runs_matrix) - \
                               special.psi(count_matrix_sequential + SMALL)  \
                               - ((-1) ** (int16(count_matrix_sequential) % 2)) / (count_matrix_sequential + 1.0) \
                               ),axis=3)

    #print "ent12_boots shape:"
    #print ent12_boots.shape
    #print "ent12_boots:      "
    #print ent12_boots
    
    mutinf_thisdof_different_sims_bootstrap_pairs = ent1_boots + ent2_boots - ent12_boots
                                                                                          
    #print "bootstrap sets:"+str(bootstrap_sets)
    #print "permutations:  "+str(permutations)
    #print "num pair runs: "+str(num_pair_runs)
    #print "mutinf sim1 sim2 shape:"

    #print mutinf_thisdof_different_sims_bootstrap_pairs.shape
    if(permutations == 0):
        print "mutinf_multinomial_difsm:"+str(average(mutinf_multinomial_sequential))
        print "stdev mutinf multi difsm:"+str(sqrt(average(var_mutinf_multinomial_sequential)))
    print "avg mutinf sim1 sim2        :"+str(average(mutinf_thisdof_different_sims_bootstrap_pairs[0,:,0]))
    mutinf_thisdof_different_sims = average(mutinf_thisdof_different_sims_bootstrap_pairs, axis=1) #average over pairs of runs in each bootstrap sample
    #now the nbins_cor*nbins_cor dimensions and num_pair_runs dimensions have been removed, leaving only bootstraps and permutations
    #print "mutinf values between sims for over original data for bootstrap=0 and permuted data\n"
    #print mutinf_thisdof_different_sims[0,0]
    #print mutinf_thisdof_different_sims[0,1:]
    if(VERBOSE >= 2):
        if bootstrap_sets > 1:
            print "mutinf values between sims for over original data for bootstrap=1 and permuted data\n"
            print mutinf_thisdof_different_sims[1,0]
            print mutinf_thisdof_different_sims[1,1:]
    #just first pair for debugging 
    #mutinf_thisdof_different_sims = mutinf_thisdof_different_sims_bootstrap_pairs[:,0]
    if(permutations > 0):
        var_ind_different_sims_pairs=zeros((bootstrap_sets,num_pair_runs),float64)
        #for bootstrap in range(bootstrap_sets):
        #   var_ind_different_sims_pairs[bootstrap,:] = vfast_var_for1Dpairs_ind(mutinf_thisdof_different_sims_bootstrap_pairs[bootstrap,:,1:])
        var_ind_different_sims_pairs = vfast_var_for1Dpairs_ind(mutinf_thisdof_different_sims_bootstrap_pairs)
        print "variance of mutinf between sims for orig data and permuts"
        print var_ind_different_sims_pairs
        var_ind_different_sims = average(var_ind_different_sims_pairs, axis=1) #average variance over pairs of runs
        #NOTE: var_ind_different_sims is of shape = (bootstrap_sets), while var_mutinf_multinomial_sequential is only 1 number
        # But this will be used in array context -- for permutations == 0, it will be implicity repeated in calc_excess_mutinf
        # for permutations > 0, it is already of shape = (bootstrap_sets)
        var_mutinf_multinomial_sequential = var_ind_different_sims 
        #print "variance of mutinf between sims for orig data and permuts"     
        
    if calc_variance == False:
       #if VERBOSE: print "   mutinf = %.5f," % (mutinf_thisdof),
       return mutinf_thisdof, Var_I , mutinf_thisdof_different_sims, var_mutinf_multinomial_sequential, average(mutinf_multinomial),average(mutinf_multinomial_sequential), pvalue, dKLtot_dKL1_dKL2, Counts_ij
    # print out the number of nonzero bins
    if VERBOSE >=2:
     #  num_nonzero_jointpop, num_nonzero_pxipyj = len(nonzero(av_joint_pop)[0]), len(nonzero(pxipyj_flat)[0]), 
     #  print "   nonzero bins (tot=%d): jointpop=%d, pxipyj=%d, combined=%d" % (nbins*nbins, num_nonzero_jointpop, num_nonzero_pxipyj, len(rnonzero_bins[nonzero_bins==True]))
       print "   mutinf this degree of freedom "+str(mutinf_thisdof)
    ##### No longer used #######################################################################
    # calculate the variance of the mutual information using its derivative and error-propagation
    #pxi_plus_pyj_flat = add.outer(chi_pop_hist1, chi_pop_hist2).flatten()
    #deriv_vector = 1 + logU[0,:] - (pxi_plus_pyj_flat[0,:]) * U
    ##############################################################################################
    
    return mutinf_thisdof, var_mi_thisdof, mutinf_thisdof_different_sims, var_mutinf_multinomial_sequential, average(mutinf_multinomial), average(mutinf_multinomial_sequential), pvalue, dKLtot_dKL1_dKL2, Counts_ij




### calculation of independent information using permutations
independent_mutinf_thisdof = None

def calc_excess_mutinf(chi_counts1, chi_counts2, bins1, bins2, chi_counts_sequential1, chi_counts_sequential2, bins1_sequential, bins2_sequential, num_sims, nbins, numangles_bootstrap,numangles, sigalpha, permutations, bootstrap_choose, calc_variance=False, which_runs=None, pair_runs=None, calc_mutinf_between_sims = "yes", file_prefix=None, plot_2d_histograms=False, adaptive_partitioning = 0):

    mutinf_tot_thisdof, var_mi_thisdof, mutinf_tot_thisdof_different_sims, var_ind_different_sims, mutinf_multinomial, mutinf_multinomial_sequential, pvalue, dKLtot_dKL1_dKL2, Counts_ij \
        = calc_mutinf_corrected(chi_counts1,chi_counts2, bins1, bins2, chi_counts_sequential1, chi_counts_sequential2, bins1_sequential, bins2_sequential, num_sims, nbins, numangles_bootstrap, numangles, calc_variance=calc_variance, bootstrap_choose=bootstrap_choose, permutations=permutations,which_runs=which_runs,pair_runs=pair_runs, calc_mutinf_between_sims=calc_mutinf_between_sims, file_prefix=file_prefix, plot_2d_histograms=plot_2d_histograms, adaptive_partitioning = adaptive_partitioning)
    
    
    
    #need to filter using p-value: for , use p-value (descriptive) to pick which terms to discard from average.
    #for bootstrap_sets=1, use p-value (Bayes) to filter mutinf values
    
    
    bootstrap_sets = mutinf_tot_thisdof.shape[0]
    sd_ind_different_sims = sqrt(var_ind_different_sims)
       
    print "mutinf_tot_thisdof shape:"+str(mutinf_tot_thisdof.shape)
    print "mutinf_tot_thisdof different sims shape:"+str(mutinf_tot_thisdof_different_sims.shape)
    print "tot_mutinfs from bootstrap samples\n:"+str(mutinf_tot_thisdof[:,0])
    if(permutations > 0 ):
        print "independent mutinf averaged over permutations\n:"+str(average(mutinf_tot_thisdof[:,1:], axis=1))
        print "independent mutinf different sims averaged over permutations\n:"+str(average(mutinf_tot_thisdof_different_sims[:,1:], axis=1))
    num_pair_runs = pair_runs.shape[0]
    independent_mutinf_thisdof                = zeros((bootstrap_sets),float64)
    independent_mutinf_thisdof_different_sims = zeros((bootstrap_sets),float64)
    if(permutations == 0):
        independent_mutinf_thisdof[:] = mutinf_multinomial
        independent_mutinf_thisdof_different_sims[:] = mutinf_multinomial_sequential
    else:
        independent_mutinf_thisdof = average(mutinf_tot_thisdof[:,1:], axis=1) # average over permutations
        independent_mutinf_thisdof_different_sims = average(mutinf_tot_thisdof_different_sims[:,1:], axis=1) #avg over permutations
    
    sd_ind_different_sims = sqrt(var_ind_different_sims)
            
    #print independent_mutinf_thisdof
    #print average(independent_mutinf_thisdof)
    #if(permutations > 0):
    #    independent_mutinf_thisdof_different_sims = average(mutinf_tot_thisdof_different_sims[:,1:],axis=1)
                
    #print "ind_mutinfs:"+str(independent_mutinf_thisdof)
    #print "tot_mutinfs_diff_sims:"+str(mutinf_tot_thisdof_different_sims)
    if(sigalpha < 1.0): #if we're doing statistics at all...
        excess_mutinf_thisdof = mutinf_tot_thisdof[:,0] - independent_mutinf_thisdof
    else:
        excess_mutinf_thisdof = mutinf_tot_thisdof[:,0]
    excess_mutinf_thisdof_different_sims = mutinf_tot_thisdof_different_sims[:,0] - independent_mutinf_thisdof_different_sims
    excess_mutinf_thisdof_different_sims -= sd_ind_different_sims
    #last term is for high pass filter, will be added back later
    #could consider having a different excess_mutinf_thisdof_different_sims for each bootstrap sample depending on the correlations between the runs it has in it
    #print "excess_mutinf_thisdof_different_sims:"+str(excess_mutinf_thisdof_different_sims)
    nonneg_excess_thisdof = logical_and((excess_mutinf_thisdof_different_sims > 0), (excess_mutinf_thisdof > 0))
    #print "nonneg excess thisdof: "+str(nonneg_excess_thisdof)
    old_excess_mutinf_thisdof = excess_mutinf_thisdof.copy()
    excess_mutinf_thisdof_different_sims += sd_ind_different_sims  #adding back the cutoff value
    excess_mutinf_thisdof[nonneg_excess_thisdof] -= excess_mutinf_thisdof_different_sims[nonneg_excess_thisdof] # subtract out high-pass-filtered mutinf for torsions in different sims
    #remember to zero out the excess mutinf thisdof from different sims that were below the cutoff value
    excess_mutinf_thisdof_different_sims[excess_mutinf_thisdof_different_sims <= sd_ind_different_sims ] = 0
    #print "corrected excess_mutinfs:"+str(excess_mutinf_thisdof)
    test_stat = 0
    mycutoff = 0
    sigtext = " "

    #now filter out those with too high of a probability for being incorrectly kept
    #print "pvalues (Bayes)"
    #print pvalue
    pvalue_toolow = pvalue > sigalpha
    if(sum(pvalue_toolow) > 0):
        print "one or more values were not significant!"
    excess_mutinf_thisdof *= (1.0 - pvalue_toolow * 1.0) #zeros elements with pvalues that are below threshold


    #dKLtot_dKL1_dKL2 *= (1.0 - pvalue_toolow * 1.0)      #zeros KLdiv Hessian matrix elements that aren't significant
    
    print "var_mi_thisdof: "+str(var_mi_thisdof)+"\n"
    print "   mutinf/ind_mutinf = cor:%.3f ex_btw:%.3f exc:%.3f ind:%.4f tot:%.3f ind_btw:%.3f tot_btw:%.3f (sd= %.3f)   %s" % (average(excess_mutinf_thisdof),  average(excess_mutinf_thisdof_different_sims), average(old_excess_mutinf_thisdof),  average(independent_mutinf_thisdof), average(mutinf_tot_thisdof[:,0]), average(independent_mutinf_thisdof_different_sims), average(mutinf_tot_thisdof_different_sims[:,0]), sqrt(average(var_mi_thisdof)),sigtext)
    #debug mutinf between sims
    #print "   mutinf/ind_mutinf = cor:%.3f ex_btw:%.3f exc:%.3f ind:%.3f tot:%.3f ind_btw:%.3f tot_btw:%.3f (sd= %.3f)   %s" % (average(excess_mutinf_thisdof),  average(excess_mutinf_thisdof_different_sims), average(old_excess_mutinf_thisdof),  average(independent_mutinf_thisdof), average(mutinf_tot_thisdof[:,0]), average(independent_mutinf_thisdof_different_sims), average(mutinf_tot_thisdof_different_sims[0,0]), var_mi_thisdof,sigtext)
    

    return excess_mutinf_thisdof, var_mi_thisdof, excess_mutinf_thisdof_different_sims, dKLtot_dKL1_dKL2, Counts_ij

# Number of chi angles per residue
NumChis = { "ALA":0, "CYS":1, "CYN":1, "CYX":1, "CY2":2, "ASP":2, "AS4":2, "ASH":2, "GLU":3, "GL4": 3, "GLH":3, "PHE":2, 
            "GLY":0, "HIS":2, "HIP":2, "HIE":2, "HID":2, "ILE":2, "LYS":4, "LYP":4, "LEU":2, "MET":3, "ASN":2, "GLN":3,
            "PRO":1, #GF: can't find def for chi2; rosetta uses chi1 only,
            "ARG":4, "SER":2, "THR":2, "VAL":1,"TRP":2, "TYR":2, "CTH":2, "F3G":1,
            "ACK":4, #acetyl-lysine   #ffAmber N and C termini next:
            "NALA":0, "NCYS":1, "NCYN":1, "NCYX":1, "NCY2":2, "NASP":2, "NAS4":2, "NASH":2, "NGLU":3, "NGL4": 3, "NGLH":3, "NPHE":2, 
            "NGLY":0, "NHIS":2, "NHIP":2, "NHIE":2, "NHID":2, "NILE":2, "NLYS":4, "NLYP":4, "NLEU":2, "NMET":3, "NASN":2, "NGLN":3,
            "NPRO":1, "NARG":4, "NSER":2, "NTHR":2, "NVAL":1,"NTRP":2, "NTYR":2,
            "CALA":0, "CCYS":1, "CCYN":1, "CCYX":1, "CCY2":2, "CASP":2, "CAS4":2, "CASH":2, "CGLU":3, "CGL4": 3, "CGLH":3, "CPHE":2, 
            "CGLY":0, "CHIS":2, "CHIP":2, "CHIE":2, "CHID":2, "CILE":2, "CLYS":4, "CLYP":4, "CLEU":2, "CMET":3, "CASN":2, "CGLN":3,
            "CPRO":1, "CARG":4, "CSER":2, "CTHR":2, "CVAL":1,"CTRP":2, "CTYR":2,
            } 


xtc_and_pdb_data = []


class XTC_wrapper:
    coords =  zeros((10, 50, 3, 200000), float32) #allocate plenty of space so it doesn't get stomped on later
    numangles = 0

xtc_coords = XTC_wrapper()

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
   dKLtot_dchis2 = []
   max_num_chis = 99
   sequential_res_num = 0
   # load angle info from xvg files
   def __str__(self): return "%s %s" % (self.name,self.num)

   def get_num_chis(self,name):
    if NumChis.has_key(name): 
        if name == "SER" or name == "THR" or name == "NSER" or name == "NTHR" or name == "CSER" or name == "CTHR":
            if(not (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(2)+name+".xvg") or os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(2)+name+".xvg.gz"))): 
                return min(self.max_num_chis, 1)
            else: 
                return min(self.max_num_chis, NumChis[name])
        else:
            return min(self.max_num_chis, NumChis[name])
    else:
        mychi = 1
        numchis = 0
        #while (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+self.num+".xvg")):
        if(not (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+".xvg") or os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+".xvg.gz"))):
            while (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+str(self.num)+".xvg") or os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+str(self.num)+".xvg.gz")  ):
                numchis += 1
            if numchis == 0:
                print "cannot find file"+self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+".xvg"
        else:
            while (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+".xvg") or os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"chi"+str(numchis+1)+name+".xvg.gz")  ):
                #quick fix until ligand mc writes residue number as well after their name
                numchis += 1
        if(numchis == 0): print "cannot find any xvg files for residue " + str(name) + str(self.num) + "\n"
        return numchis  

   def has_phipsi(self,name):
     has_phipsi = False
     if NumChis.has_key(name): 
         has_phipsi = True
     else:
         if (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"phi"+name+".xvg")):
            has_phipsi = True
            if (os.path.exists(self.xvg_basedir+"run1"+self.xvg_chidir+"psi"+name+".xvg")):
                has_phipsi = True
            else:
                has_phipsi = False
         else:
            has_phipsi = False
     return has_phipsi

   def _load_xvg_data(self, basedir, num_sims, max_angles, chi_dir = "/dihedrals/g_chi/", skip=1, skip_over_steps=0,coarse_discretize=None):
      myname = self.name
      mynumchis = self.get_num_chis(myname)
      shifted_angles = zeros((self.nchi,num_sims,max_angles),float64)
      #self.numangles[:] = run_params.num_structs
            
      #shifted_angles[:,:,:] = -999 #a null value other than zero
      #weird fix for residue type "CYS2" in Gromacs
      #if myname == "CYS": myname += "2"
      #assert num_sims == len(self.which_runs[0])
          
      res_num = str(self.xvg_resnum)

      if(coarse_discretize == None ):
       for chi_num in range(self.nchi):
         #print "Chi:"+str(chi_num+1)+"\n"
         for sequential_sim_num in range(num_sims):
            #sim_index_str = str(self.which_runs[sequential_sim_num])
            if(self.nchi == mynumchis + 2): #phi/psi
                if(chi_num == 0):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname+res_num+".xvg"
                if(chi_num == 1):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname+res_num+".xvg"
                if(chi_num > 1 and chi_num <= mynumchis + 1):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname+res_num+".xvg"
            else:
                if(chi_num < mynumchis):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname+res_num+".xvg"
                if(chi_num == mynumchis):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname+res_num+".xvg"
                if(chi_num == mynumchis + 1):
                    xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname+res_num+".xvg"
            # gf: gromacs is changing some of my residue names when I use it to read in a PDB trajectory
            if not (os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz")):
               if myname == "LYS":
                  for myname_alt in ("LYSH", "LYP"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break
               if myname == "ASP" or myname == "ASH" or myname == "AS4":
                  for myname_alt in ("ASH","AS4", "ASP"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break   
               if myname == "GLU" or myname == "GLH" or myname == "GL4":
                  for myname_alt in ("GLH","GL4","GLU"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break   
               if myname == "HIS" or myname == "HID" or myname == "HIE" or myname == "HID" or myname == "HIP":
                  for myname_alt in ("HISA", "HISB", "HIE", "HID", "HIS", "HIP"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break
	       if myname == "CYS":
		  for myname_alt in ("CYS2", "CYX", "CYN", "F3G"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break
               if myname == "ALA":
		  for myname_alt in ("NALA", "NAL"):
                      if(self.nchi == mynumchis + 2): #phi/psi
                          if(chi_num == 0):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                          if(chi_num > 1 and chi_num <= mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num-1)+myname_alt+res_num+".xvg"
                      else:
                          if(chi_num < mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"chi"+str(chi_num+1)+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg"
                          if(chi_num == mynumchis + 1):
                              xvg_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"
                      if os.path.exists(xvg_fn) or os.path.exists(xvg_fn+".gz"): break

            if os.path.exists(xvg_fn):
                (data,titlestr)=readxvg(xvg_fn,skip,skip_over_steps)
            elif os.path.exists(xvg_fn+".gz"):
                (data,titlestr)=readxvg(xvg_fn+".gz",skip,skip_over_steps)
            else:
                print "ERROR: unable to find file '%s[.gz]'" % xvg_fn
                sys.exit(1)
                
            self.numangles[sequential_sim_num] = len(data[:,1])
            targ_shape = shape(self.angles)
            data_shape = len(data[:,1])
            #print "sim_num "+str(sequential_sim_num)+" numangles: "+str(self.numangles[sequential_sim_num])+" shape of target array: "+str(targ_shape)+" shape of data array: "+str(data_shape)
            data[:,1] = (data[:,1] + 180)%360 - 180 #Check and make sure within -180 to 180; I think this should do it 
            self.angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = data[:,1]
      
      else: #coarse_discretize 
             for sequential_sim_num in range(num_sims):
                    # gf: gromacs is changing some of my residue names when I use it to read in a PDB trajectory
                    xvg_fn_phi = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname+res_num+".xvg"
                    xvg_fn_psi = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname+res_num+".xvg"
                    def xvg_filenames_phipsi(thisname):
                           return [basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname_alt+res_num+".xvg", \
                                   basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname_alt+res_num+".xvg"]
                                         
                    if not (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                           (os.path.exists(xvg_fn_psi) or os.path.exists(xvg_fn_psi+".gz")) :
                           if myname == "LYS":
                                  for myname_alt in ("LYSH", "LYP"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt) 
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                           if myname == "ASP" or myname == "ASH" or myname == "AS4":
                                  for myname_alt in ("ASH","AS4", "ASP"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt)
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                           if myname == "GLU" or myname == "GLH" or myname == "GL4":
                                  for myname_alt in ("GLH","GL4","GLU"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt)
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                           if myname == "HIS" or myname == "HID" or myname == "HIE" or myname == "HID" or myname == "HIP":
                                  for myname_alt in ("HISA", "HISB", "HIE", "HID", "HIS", "HIP"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt)
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                           if myname == "CYS":
                                  for myname_alt in ("CYS2", "CYX", "CYN", "F3G"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt)
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                           if myname == "ALA":
                                  for myname_alt in ("NALA", "NAL"):
                                         (xvg_fn_phi, xvg_fn_psi) = xvg_filenames_phipsi(myname_alt)
                                         if (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")) and \
                                            (os.path.exists(xvg_fn_phi) or os.path.exists(xvg_fn_phi+".gz")): break
                    else:
                           xvg_fn_phi = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"phi"+myname+res_num+".xvg"
                           xvg_fn_psi = basedir+"run"+str(sequential_sim_num+1)+chi_dir+"psi"+myname+res_num+".xvg"
                    print "pathname phi: "
                    print xvg_fn_phi
                    print "pathname psi: "
                    print xvg_fn_psi
                    if os.path.exists(xvg_fn_phi):
                           (data_phi,titlestr)=readxvg(xvg_fn_phi,skip,skip_over_steps)
                    elif os.path.exists(xvg_fn_phi+".gz"):
                           (data_phi,titlestr)=readxvg(xvg_fn_phi+".gz",skip,skip_over_steps)
                    else:
                           print "ERROR: unable to find file '%s[.gz]'" % xvg_fn_phi
                           sys.exit(1)
                    if os.path.exists(xvg_fn_psi):
                           (data_psi,titlestr)=readxvg(xvg_fn_psi,skip,skip_over_steps)
                    elif os.path.exists(xvg_fn_phi+".gz"):
                           (data_psi,titlestr)=readxvg(xvg_fn_psi+".gz",skip,skip_over_steps)
                    else:
                           print "ERROR: unable to find file '%s[.gz]'" % xvg_fn_psi
                           sys.exit(1)
                    self.numangles[sequential_sim_num] = len(data_phi[:,1])   
                    data_phi[:,1] = (data_phi[:,1] + 180)%360 - 180 #Check and make sure within -180 to 180; I think this should do it
                    data_psi[:,1] = (data_phi[:,1] + 180)%360 - 180 #Check and make sure within -180 to 180; I think this should do it
                    #print data_phi
                    for i in range(self.numangles[sequential_sim_num]):
                          phi = data_phi[i,1]
                          psi = data_psi[i,1]
                          ## Alpha   ... these angle numbers are just fixed values...
                          if ( -180 < phi < 0 and -100 < psi < 45):
                                 self.angles[0,sequential_sim_num,i] = -179.9
                                 #print "alpha"
                          ## Beta
                          elif ( -180 < phi < -45 and (45 < psi or psi > -135) ):
                                 self.angles[0,sequential_sim_num,i] = -89.9
                                 #print "beta"
                          ## Turn
                          elif ( 0 < phi < 180 and -90 < psi < 90):
                                 self.angles[0,sequential_sim_num,i] = 0.1
                                 #print "turn"
                          ## Other
                          else:  
                                 self.angles[0,sequential_sim_num,i] = 90.1
                                 #print "other"
                    #print "self.angles:"+str(self.angles[0,sequential_sim_num])
                    #print self.angles[0,sequential_sim_num]

          
      for sequential_sim_num in range(self.num_sims): ## PATCH to fix problem with different # of angles in different sims
              self.numangles[sequential_sim_num] = min(self.numangles[:])     ## PATCH to fix problem with different # of angles in different sims

# for phi/psi discretization, use made-up numbers like 45, 135, 215, etc. and one torsion per residue, even making temporary file if needed
##


   
  
   def load_xtc_data(self, basedir, num_sims, max_angles, chi_dir = "/dihedrals/g_chi/", skip=1, skip_over_steps=0, pdbfile=None, xtcfile=None):
          
          #dataarray=zeros((int((inlength-skip_over_steps)/skip) + extra_record,numfields),float64) #could start with zero, so add + extra_record
          
      skiplines=0
   #Read data into array
      #for i in range(int((inlength-skip_over_steps)/skip) + extra_record): #could start with zero, so add + extra_record ...
      #   if(i*skip + skip_over_steps < inlength): # ... but make sure we don't overshoot
      #       entries=inlines[i*skip+skip_over_steps].split()
      global xtc_and_pdb_data
      global xtc_coords
      myname = self.name
      self.nchi = 3
      mynumchis = 3 #self.get_num_chis(myname)
      tot_residues = 0
      #shifted_angles = zeros((self.nchi,num_sims,max_angles),float64)
      #self.numangles[:] = run_params.num_structs
      #shifted_angles[:,:,:] = -999 #a null value other than zero
      #weird fix for residue type "CYS2" in Gromacs
      #if myname == "CYS": myname += "2"
      #assert num_sims == len(self.which_runs[0])
      #res_num = str(self.sequential_res_num + 1) #since MDAnalysis starts with atom 1 
      if(skip == 1):
             extra_record = 0
      else:
             extra_record = 1
      print "sequential residue number: "+str(self.sequential_res_num)
      for sequential_sim_num in range(num_sims):
             xtc_fn = basedir+"run"+str(sequential_sim_num+1)+chi_dir+xtcfile+str(sequential_sim_num)+".xtc"
             pdb_fn = pdbfile
             if os.path.exists(xtc_fn) and os.path.exists(pdb_fn):
                 if (xtc_and_pdb_data == []):
                     i = 0
                     tempdata = MDAnalysis.Universe(pdb_fn, xtc_fn)
                     tot_residues = shape(tempdata.atoms.coordinates())[0]
                     print "total residues:"+str(tot_residues)
                     count = 0
                     for ts in tempdata.trajectory:
                            if(count > skip_over_steps): 
                                   if( count % skip == 0):
                                          i += 1
                            count += 1
                     xtc_coords.numangles = i
                     self.numangles[sequential_sim_num] = i
                     xtc_coords.coords = resize(xtc_coords.coords,(num_sims, tot_residues, 3, self.numangles[0])) #shrink down to proper limits
                     tempdata.trajectory.close_trajectory()
                 if len(xtc_and_pdb_data) < num_sims: #if we haven't opened all the trajectory files yet
                     xtc_and_pdb_data.append(MDAnalysis.Universe(pdb_fn, xtc_fn))
                     if(len(xtc_and_pdb_data) >= 1):
                         i = 0
                         count = 0
                         print "loading cartesian data, run: "+str(sequential_sim_num)
                         for ts in xtc_and_pdb_data[sequential_sim_num].trajectory:
                                #print shape((xtc_and_pdb_data[sequential_sim_num].atoms.coordinates())[:,:])
                                #print shape(xtc_coords.coords[sequential_sim_num, :tot_residues, :, i])
                                if(count > skip_over_steps): 
                                   if( count % skip == 0):
                                          xtc_coords.coords[sequential_sim_num, :tot_residues, :3, i] = (xtc_and_pdb_data[sequential_sim_num].atoms.coordinates())[:tot_residues,:3]
                                          i += 1
                                count += 1
                         self.numangles[sequential_sim_num] = i
                         xtc_and_pdb_data[sequential_sim_num].trajectory.rewind() #reset for next use
                 else: #get numangles from xtc_coords
                     self.numangles[sequential_sim_num] = xtc_coords.numangles
                 #print shape(xtc_coords.coords)
                 #print xtc_coords.coords[sequential_sim_num, self.sequential_res_num, :, :self.numangles[sequential_sim_num]]
                 #print xtc_coords.coords[sequential_sim_num, :, 0, :self.numangles[sequential_sim_num]]
                 #print xtc_coords.coords[0, :, 0, :]
                 
                 
                 
                 
                 self.angles[:3,sequential_sim_num,:self.numangles[sequential_sim_num]] =  xtc_coords.coords[sequential_sim_num, self.sequential_res_num, :, :self.numangles[sequential_sim_num]]
                 #if( i % 100 == 0): 
                 #print "timestep "+str(i)
                 
                 
             else:
                 if os.path.exists(xtc_fn):
                     print "ERROR: unable to find file '%s[.gz]'" % pdb_fn
                 else:
                     print "ERROR: unable to find file '%s[.gz]'" % xtc_fn
                 sys.exit(1)
      for sequential_sim_num in range(self.num_sims): ## PATCH to fix problem with different # of angles in different sims
             self.numangles[sequential_sim_num] = min(self.numangles[:])     ## PATCH to fix problem with different # of angles in different sims
      
      ## RESCALE CARTESIANS 
      
      ## first, zero centroid of all particles
      ## ideally filter out the first eigenvector, assuming rot+trans alignment already performed, so this is in effect filters out #7
      
      ##Get coordinate ranges for all 3 cartesians over all residues
      #Perhaps should have box sizes grabbed somehow from xtc file? 
      #coordmin = amin(amin(self.angles[:3,:,:min(self.numangles)],axis=2),axis=1)
      #coordrange = amax(amax(self.angles[:3,:,:min(self.numangles)],axis=2),axis=1) - coordmin + SMALL #+ SMALL to avoid div by zero
      #coordavg = average(average(self.angles[:3,:,:min(self.numangles)],axis=2),axis=1)
   
      # perform rescaling so that cartesians are in range (-120, 120), this will give wiggle room for substantial shape changes
      #for sequential_sim_num in range(self.num_sims):
      #    offset = zeros((3,self.numangles[sequential_sim_num]), float64)
      #    scaling_factor = ones((3,self.numangles[sequential_sim_num]), float64)
      #    for mycoord in range(3):
      #        offset[mycoord,:] = coordavg[mycoord]
      #        scaling_factor = (240 - SMALL) / (coordrange[mycoord])
      #    self.angles[:3,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
      #        (self.angles[:3,sequential_sim_num,:self.numangles[sequential_sim_num]]  -  offset) * scaling_factor - 120 
      
             

                 
          



       


      ### Rank-order data for adaptive partitioning -- note that I give the rank over the pooled data for each angle in each sim
##      shifted_angles = self.angles.copy()
##     shifted_angles[shifted_angles < 0] += 360 # wrap these around so that -180 and 180 are part of the same bin
##      for chi_num in range(self.nchi):
##         #print "Chi:"+str(mychi+1)+"\n"
##         shifted_flat = resize(((swapaxes((shifted_angles[chi_num,:,:]).copy(),0,1))),(sum(self.numangles)))
##         sorted_angles = sort(shifted_flat)
##         #print "numangles: "+str(self.numangles)+" sum_numangles: "+str(sum(self.numangles))+" num of sorted angles: "+str(len(sorted_angles))+" num of datapoints: "+str(len(shifted_flat))
##         #print (searchsorted(sorted_angles,shifted_angles[chi_num,0,:]))
##         for sequential_sim_num in range(num_sims):
##             #print "Shifted Angles: Length: "+str(shape(shifted_angles[chi_num,0,:]))+"\n"
##             #print shifted_angles[chi_num,sim_num,:]
##             #print "Sorted Angles: Length\n"+str(shape(sorted_angles))+"\n"
##             #print sorted_angles[:]
##             #print "Rank ordered Angles \n"
##            self.rank_order_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
##                                      searchsorted(sorted_angles,shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles
             #print self.rank_order_angles[chi_num,sim_num,:]
	     #print "\n"
   
   # load angle info from an all_angle_info object
   def _load_pdb_data(self, all_angle_info, max_angles):
      #shifted_angles = zeros((self.nchi,all_angle_info.num_sims,max_angles),float64)
      #shifted_angles[:,:,:] = -999 #a null value other than zero      
      for sequential_sim_num in range(self.num_sims): #range(all_angle_info.num_sims):
          curr_angles, numangles = all_angle_info.get_angles(sequential_sim_num, int(self.xvg_resnum), self.name, self.backbone_only, self.phipsi, self.max_num_chis)
          self.angles[:, sequential_sim_num, 0:numangles] = (curr_angles + 180)%360 - 180
          self.numangles[sequential_sim_num]= numangles
      for sequential_sim_num in range(self.num_sims): ## PATCH to fix problem with different # of angles in different sims
          self.numangles[sequential_sim_num] = min(self.numangles[:])     ## PATCH to fix problem with different # of angles in different sims
      ### Rank-order data for adaptive partitioning -- note that I give the rank over the pooled data for each angle in each sim



##      shifted_angles = self.angles.copy()
##      shifted_angles[shifted_angles < -90] += 360 # wrap these around so that -180 and 180 are part of the same bin
##      for chi_num in range(self.nchi):
##         #print "Chi:"+str(mychi+1)+"\n"
##         shifted_flat = resize(((swapaxes((shifted_angles[chi_num,:,:]).copy(),0,1))),(sum(self.numangles)))
##         sorted_angles = sort(shifted_flat)
##         #print "numangles: "+str(self.numangles)+" sum_numangles: "+str(sum(self.numangles))+" num of sorted angles: "+str(len(sorted_angles))+" num of datapoints: "+str(len(shifted_flat))
##         #print shifted_angles[chi_num,0,:200]
##         #print sorted_angles[:200]
##         #print (searchsorted(sorted_angles,shifted_angles[chi_num,0,:]))
##         for sequential_sim_num in range(self.num_sims):
##             self.rank_order_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
##                                   searchsorted(sorted_angles,shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles
             #print self.rank_order_angles[chi_num,sim_num,:]
            
   # load the angles and calculate the total entropy of this residue and its variance.
   # The entropy of a residue with two chi angles is calculated as
   #     H(X1,X2) = H(X1) + H(X2) - I(X1,X2)
   # where I() is the mutual information between X1 and X2.
   # For residues with more than 3 chi angles, we assume the 2nd order approximation is sufficient.
   # The variance residue's entropy is the sum of the variances of all the terms.

   def correct_and_shift_carts(self,num_sims,bootstrap_sets,bootstrap_choose):
       shifted_carts = self.angles.copy()
       myname = self.name
       mynumchis = self.get_num_chis(myname)
       # shift cartesians to [0, 360] so binning will be compatible for angular data
       # we don't know if all the sims are aligned togther.. if so, then shifting shoud be done for whole dataset not each sim separately 
       for i in range((shape(shifted_carts))[0]):
              for j in range ((shape(shifted_carts))[1]):
                     shifted_carts[i,j,:] = shifted_carts[i,j,:] - average(shifted_carts[i,j,:])
                     shifted_carts[i,j,:] = shifted_carts[i,j,:] * 360.0 / \
                                            (max(shifted_carts[i,j,:]) - min(shifted_carts[i,j,:]))
       
       self.angles[:,:,:] = shifted_carts[:,:,:]
       #del shifted_carts

   def expand_contract_data(self,num_sims,bootstrap_sets,bootstrap_choose):
       shifted_data = self.angles.copy()
       myname = self.name
       mynumchis = self.get_num_chis(myname)
       mymin = self.minmax[0,:]
       mymax = self.minmax[1,:]
       mymin_array1 = zeros((num_sims), float64)
       mymin_array2 = zeros((num_sims, shape(self.angles)[-1]), float64)
       mymax_array1 = zeros((num_sims), float64)
       mymax_array2 = zeros((num_sims, shape(self.angles)[-1]), float64)
       
       # shift data to [0, 360] so binning will be compatible for max and min of data
       # have center of data be max - min / 2 rather than average
       for mychi in range(self.nchi):
              mymin_array1[:] = mymin[mychi]
              mymin_array2[:, :] = resize(mymin_array1,shape(mymin_array2))
              mymax_array1[:] = mymax[mychi]
              mymax_array2[:, :] = resize(mymax_array1,shape(mymax_array2))

              #midpoint = ( mymax_array2[:,:] - mymin_array2[:,:] ) / 2.0 
              zeropoint =  mymin_array2[:,:] 
              # no recentering, as I think it might introduce bias
              # rather, we just stretch the range by setting the min to zero,
              # then grabbing on the max and pulling it to 360
              shifted_data[mychi,:,:] = shifted_data[mychi,:,:] - zeropoint 
              shifted_data[mychi,:,:] = shifted_data[mychi,:,:] * 360.0 / \
                                            ( mymax_array2[:,:] - mymin_array2[:,:] ) #+ midpoint ##if using midpoint recentering
              shifted_data[mychi,:,:] = (shifted_data[mychi,:,:])%360  #Check and make sure within 0 to 360
            
       for mysim in range(num_sims):
              shifted_data[:,mysim,self.numangles[mysim]:] = 0
       print "orig data:"
       print self.angles[mychi,0,:min(self.numangles)]
       print "data range:"
       print mymin_array1
       print mymax_array1
       print "reshaped data:"
       print shifted_data[mychi,0,:min(self.numangles)]
       self.angles[:,:,:] = shifted_data[:,:,:]
       #del shifted_data
   
   def correct_and_shift_angles(self,num_sims,bootstrap_sets,bootstrap_choose, coarse_discretize = None):
       ### Rank-order data for adaptive partitioning --
       ### note that I give the rank over the pooled data for angles from all sims
   
       myname = self.name
       mynumchis = self.get_num_chis(myname)
       print "residue name: "+str(self.name)+" num chis: "+str(mynumchis)+"\n"
       shifted_angles = self.angles.copy()
       shifted_angles[shifted_angles < 0] += 360 # wrap these around so that -180 and 180 are part of the same bin
       assert(self.nchi > 0)
       
       #wrap torsions of 2-fold symmetric and 3-fold symmetric terminal chi angles
       # ARG and LYS's "chi5" is symmetric, but these isn't one of the standard chi angles
       # and typically we only have up to chi4 for these anyways.

       #However, we do not correct protonated ASP or GLU residues, as binding of a proton to these breaks symmetry.

       if(CORRECT_FOR_SYMMETRY == 1):
        if(self.nchi == mynumchis + 2): #phi/psi            
               if (myname == "ASP" or myname == "GLU" or myname == "PHE" or myname == "TYR" \
                   or myname == "NASP" or myname == "NGLU" or myname == "NPHE" or myname == "NTYR" \
                   or myname == "CASP" or myname == "CGLU" or myname == "CPHE" or myname == "CTYR") and coarse_discretize == None:
                      #last chi angle is 2-fold symmetric
                      myangles = shifted_angles[mynumchis + 1,:,:]
                      myangles[myangles > 180] = myangles[myangles > 180] - 180
                      shifted_angles[mynumchis + 1,:,:] = myangles
                      self.symmetry[mynumchis + 1] = 2
              
        else:
            if(self.nchi == mynumchis):
                   if (myname == "ASP" or myname == "GLU" or myname == "PHE" or myname == "TYR" \
                      or myname == "NASP" or myname == "NGLU" or myname == "NPHE" or myname == "NTYR" \
                      or myname == "CASP" or myname == "CGLU" or myname == "CPHE" or myname == "CTYR") and coarse_discretize == None :
                          #last chi angle is 2-fold symmetric
                          myangles = shifted_angles[mynumchis - 1,:,:]
                          myangles[myangles > 180] = myangles[myangles > 180] - 180
                          shifted_angles[mynumchis - 1,:,:] = myangles
                          self.symmetry[mynumchis - 1] = 2
            
              
       self.angles[:,:,:] = shifted_angles[:,:,:]  #now actually shift the angles, important for entropy and non-adaptive partitioning correlations
       
       #unshift 
       #self.angles[self.angles > 180] -= 360 # back to [-180, 180]
       for chi_num in range(self.nchi):
         #print "Chi:"+str(mychi+1)+"\n"
         shifted_flat = resize(((swapaxes((shifted_angles[chi_num,:,:]).copy(),0,1))),(sum(self.numangles)))
         self.sorted_angles[chi_num,:] = sort(shifted_flat)
         #print "numangles: "+str(self.numangles)+" sum_numangles: "+str(sum(self.numangles))+" num of sorted angles: "+str(len(sorted_angles))+" num of datapoints: "+str(len(shifted_flat))
         #print (searchsorted(sorted_angles,shifted_angles[chi_num,0,:]))
         for sequential_sim_num in range(num_sims):
             #print "Shifted Angles: Length: "+str(shape(shifted_angles[chi_num,0,:]))+"\n"
             #print shifted_angles[chi_num,sim_num,:]
             #print "Sorted Angles: Length\n"+str(shape(sorted_angles))+"\n"
             #print sorted_angles[:]
             #print "Rank ordered Angles \n"
             self.rank_order_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
              searchsorted(self.sorted_angles[chi_num,:], shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles for all sims together
             sorted_angles_sequential = sort((shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]).copy())
             self.rank_order_angles_sequential[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]] = \
                                      searchsorted(sorted_angles_sequential, shifted_angles[chi_num,sequential_sim_num,:self.numangles[sequential_sim_num]]) # rank-ordered dihedral angles, each sim ranked separately ... for corr between sims                                                              
             #print self.rank_order_angles[chi_num,sim_num,:]
	     #print "\n"
         for bootstrap in range(bootstrap_sets):
             boot_numangles  = 0
             for boot_sim_num in range(bootstrap_choose):
                 sequential_sim_num = self.which_runs[bootstrap,boot_sim_num]
                 boot_numangles += self.numangles[sequential_sim_num]
             boot_angles = zeros((bootstrap_choose,boot_numangles),float64)
             for boot_sim_num in range(bootstrap_choose):
                 #sequential_sim_num = self.which_runs[bootstrap,boot_sim_num]
                 #copy angles
                 #self.boot_sorted_angles[chi_num,bootstrap,:]
                 boot_angles[boot_sim_num,:self.numangles[sequential_sim_num]]= self.angles[chi_num, self.which_runs[bootstrap,boot_sim_num], :self.numangles[sequential_sim_num]]
             # resize for sorting so the extra angle slots go into axis 0 and so will appear after all the data
             boot_flat = resize(swapaxes(boot_angles,0,1),(boot_numangles))
             #print shape(boot_flat)
             self.boot_sorted_angles[chi_num,bootstrap,:boot_numangles] = sort(boot_flat)
             self.boot_ranked_angles[chi_num,bootstrap,:boot_numangles] = searchsorted( \
                 self.boot_sorted_angles[chi_num,bootstrap,:boot_numangles], \
                 boot_flat)
             self.numangles_bootstrap[bootstrap] = boot_numangles
             
       return

   
   def __init__(self,myname,mynum,xvg_resnum,basedir,num_sims,max_angles,xvgorpdb,binwidth,sigalpha=1,
                permutations=0,phipsi=0,backbone_only=0,adaptive_partitioning=0,which_runs=None,pair_runs=None,bootstrap_choose=3,
                calc_variance=False, all_angle_info=None, xvg_chidir = "/dihedrals/g_chi/", skip=1, skip_over_steps=0, calc_mutinf_between_sims="yes", max_num_chis=99,
                sequential_res_num = 0, pdbfile = None, xtcfile = None, output_timeseries = "no", minmax = None, bailout_early = False):
      global xtc_coords 
      self.name = myname
      self.num = mynum
      self.xvg_basedir = basedir
      self.xvg_chidir = xvg_chidir
      self.xvg_resnum = xvg_resnum
      self.sequential_res_num = sequential_res_num
      self.backbone_only, self.phipsi = backbone_only, phipsi
      self.max_num_chis = max_num_chis
      coarse_discretize = None
      if(phipsi >= 0): 
             self.nchi = self.get_num_chis(myname) * (1 - backbone_only) + phipsi * self.has_phipsi(myname)
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
      check_for_free_mem()
      
      #allocate stuff
      bootstrap_sets = self.which_runs.shape[0]
      self.entropy =  zeros((bootstrap_sets,self.nchi), float64)
      self.entropy2 =  zeros((bootstrap_sets,self.nchi), float64) #entropy w/fewer bins
      self.entropy3 =  zeros((bootstrap_sets,self.nchi), float64) #entropy w/fewer bins
      self.var_ent =  zeros((bootstrap_sets,self.nchi), float64)
      self.numangles_bootstrap = zeros((bootstrap_sets),int32)
      print "\n#### Residue: "+self.name+" "+self.num+" torsions: "+str(self.nchi), utils.flush()
      binwidth = float(binwidth)
      bins = arange(0,360, binwidth) #  bin edges global variable
      nbins=len(bins) # number of bins
      nbins_cor = int(nbins * FEWER_COR_BTW_BINS);
      sqrt_num_sims=sqrt(num_sims)
      self.chi_pop_hist=zeros((bootstrap_sets, self.nchi,nbins),float64)
      self.chi_counts=zeros((bootstrap_sets, self.nchi, nbins), int32)
      #self.chi_var_pop=zeros((bootstrap_sets, self.nchi,nbins),float64)
      self.chi_pop_hist_sequential=zeros((num_sims, self.nchi, nbins_cor), float64)
      num_histogram_sizes_to_try = 2  # we could try more and pick the optimal size
      self.chi_counts_sequential=zeros((num_sims, self.nchi, nbins_cor), int32) #half bin size
      self.chi_counts_sequential_varying_bin_size=zeros((num_histogram_sizes_to_try, num_sims, self.nchi, int(nbins*(num_histogram_sizes_to_try/2)) ), int32) #varying bin size
      self.angles = zeros((self.nchi,num_sims,max_angles),float64)         # the dihedral angles
      #self.sorted_angles = zeros((self.nchi,num_sims,max_angles),float64) # the dihedral angles sorted
      self.boot_sorted_angles = zeros((self.nchi,bootstrap_sets,bootstrap_choose*max_angles),float64)
      self.boot_ranked_angles = zeros((self.nchi,bootstrap_sets,bootstrap_choose*max_angles),int32) 
      self.rank_order_angles = zeros((self.nchi,num_sims,max_angles),int32) # the dihedral angles
                                                                        # rank-ordered with respect to all sims together
      self.rank_order_angles_sequential = zeros((self.nchi,num_sims,max_angles),int32) # the dihedral angles
                                                                        # rank-ordered for each sim separately
                                                                        # for mutinf between sims
      #self.ent_hist_left_breaks = zeros((self.nchi, nbins * MULT_1D_BINS + 1),float64)
      #self.ent_hist_binwidths = zeros((bootstrap_sets, self.nchi, nbins * MULT_1D_BINS),float64)
      self.ent_from_sum_log_nn_dists = zeros((bootstrap_sets, self.nchi, MAX_NEAREST_NEIGHBORS),float64)
      self.minmax = zeros((2,self.nchi))
      self.minmax[1,:] += 1 #to avoid zero in divide in expand_contract angles
      ### load DATA
      if(xvgorpdb == "xvg"):
         self._load_xvg_data(basedir, num_sims, max_angles, xvg_chidir, skip,skip_over_steps, coarse_discretize)
      if(xvgorpdb == "pdb"):
         self._load_pdb_data(all_angle_info, max_angles)
      if(xvgorpdb == "xtc"):
         self.load_xtc_data(basedir, num_sims, max_angles, xvg_chidir, skip, skip_over_steps, pdbfile, xtcfile) 

      if ((self.name in ("GLY", "ALA")) and phipsi == 0): return

      self.sorted_angles = zeros((self.nchi, sum(self.numangles)),float64)

      if(xvgorpdb == "xvg" or xvgorpdb == "pdb"):
             self.correct_and_shift_angles(num_sims,bootstrap_sets,bootstrap_choose, coarse_discretize)
      
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
      
      if(xvgorpdb == "xtc"):
             self.correct_and_shift_carts(num_sims,bootstrap_sets,bootstrap_choose)
      
      if(bailout_early == True): #if we just want the angles and especially the min and max values without the binning, etc....
             print "bailing out early with min/max values"
             return
      

      inv_binwidth = 1.0 / binwidth
      inv_binwidth_cor = nbins_cor / 360.0;
      #print "binwidth:"+str(binwidth)
      #if(adaptive_partitioning == 1):
      #print "numangles:"+str((self.numangles))+"sum numangles: "+str(sum(self.numangles))+" binwidth = "+str(sum(self.numangles)/nbins)+"\n"
      inv_binwidth_adaptive_bootstraps = (nbins * 1.0 /self.numangles_bootstrap)  #returns a vector
      inv_binwidth_adaptive_sequential = (nbins_cor * 1.0 /(self.numangles))  #returns a vector
      #print "inv_binwidth_adaptive: "+str(inv_binwidth_adaptive)

      number_per_ent_bin = sum(self.numangles) / (nbins * MULT_1D_BINS)

              
      #print shape(self.ent_hist_left_breaks)
      #print shape(self.ent_hist_binwidths)
      #self.ent_pdf_adaptive=zeros((bootstrap_sets,self.nchi,MULT_1D_BINS * nbins),float64) # normalized number of counts per bin,
      #for bootstrap in range(bootstrap_sets):
      #for mychi in range(self.nchi):
      #    for i in range(nbins*MULT_1D_BINS):   #use sorted angles. Note, we've only sorted angles once, could be done for each bootstrap sample just for the 1D entropies -- sorting each bootstrap sample separately would confuse correlations in rank order space
      #        self.ent_hist_left_breaks[mychi,i] = self.sorted_angles[mychi, number_per_ent_bin * i]
      #    self.ent_hist_left_breaks[mychi,(nbins*MULT_1D_BINS + 1) - 1] = self.sorted_angles[mychi,sum(self.numangles) - 1] + 0.0001
      #    for i in range(nbins*MULT_1D_BINS):   #use sorted angles. Note, we've only sorted angles once, could be done for each bootstrap sample just for the 1D entropies -- sorting each bootstrap sample separately would confuse correlations in rank order space
                  #if(i == nbins - 1):
      #            self.ent_hist_binwidths[:,mychi, i]  = resize((TWOPI / 360.0) * (self.sorted_angles[mychi, -1] - self.ent_hist_left_breaks[mychi,i]) / self.symmetry[mychi], self.ent_hist_binwidths.shape[0])
                  #else:
      #            self.ent_hist_binwidths[:,mychi, i]  = resize((PI / 180.0) * (self.ent_hist_left_breaks[mychi, i + 1] - self.ent_hist_left_breaks[mychi, i]), self.ent_hist_binwidths.shape[0])




      




      #print "Average Number of Angles Per Entropy Bin:\n"
      #print number_per_ent_bin
      #print "Entropy Histogram Left Breaks:\n"
      #print self.ent_hist_left_breaks[:,:]
      #print "Entropy Histogram Binwidths:\n"
      #print self.ent_hist_binwidths[0,:,:]
      #print "Sum of Histogram Binwidhths:\n"
      #print sum(self.ent_hist_binwidths,axis=-1)
      #use these left histogram breaks and binwidths below
      
      if VERBOSE >= 2: print "Angles: ", map(int, list(self.angles[0,0,0:self.numangles[0]])), "\n\n"
      
      
      # find the bin for each angle and the number of counts per bin
      # need to weave this loop for speed
      max_num_angles = int(max(self.numangles))
      self.bins = zeros((self.nchi, self.permutations + 1, bootstrap_sets, bootstrap_choose * max_num_angles ), int8) # the bin for each dihedral
      self.simbins = zeros((self.nchi, self.permutations + 1, num_sims, max_num_angles), int8)
      self.counts=zeros((bootstrap_sets,self.nchi,MULT_1D_BINS * nbins),int32) # number of counts per bin
      self.counts2=zeros((bootstrap_sets,self.nchi, nbins),int32) # number of counts per bin w/fewer bins
      self.counts3=zeros((bootstrap_sets,self.nchi, SAFE_1D_BINS),int32) # number of counts per bin w/ even fewer bins
      
      #counts_marginal=zeros((bootstrap_sets,self.nchi,nbins),float32) # normalized number of counts per bin, 
      counts_adaptive=zeros((bootstrap_sets,self.nchi,MULT_1D_BINS * nbins),int32)
      
#
#
#work zone for weaved replacement of loops below
#
#
#def binsingle(angle,inv_binwidth):
#   if angle == -180: angle = 180
#   return int(floor((angle-0.00001 + 180)*inv_binwidth)) #so we don't get an overshoot if angle is exactly 180

#def binsingle_adaptive(angle,inv_binwidth):
    #print "rank: "+str(angle)+" binwidth: "+str(1.0/inv_binwidth)+" bin: "+str(int(floor(angle*inv_binwidth)))
#    return int(floor(angle*inv_binwidth)) #here "angle" is a rank-order for the angle over sum(numangles)


#                    printf("simnum %d mynumangles %d\\n",simnum,mynumangles);
#                                          printf("bin %d \\n",bin1);
      code_nonadaptive_doublebins = """
              // weave8
              int mynumangles, mynumangles_sum, bin1, bin2, bin3, simnum; //bin2 and bin3 are for lower-res histograms
              double angle;
              for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                mynumangles = 0;
                *(numangles_bootstrap + mybootstrap) = 0 ;   
                for(int mysim=0; mysim < bootstrap_choose; mysim++) {
                    simnum = *(which_runs + mybootstrap*bootstrap_choose + mysim);
                    mynumangles = *(numangles + simnum);
                    mynumangles_sum = *(numangles_bootstrap + mybootstrap);

                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(angles + mychi*num_sims*max_angles  + simnum*max_angles + anglenum);
                      if(angle > 360) angle = angle - 360;
                      if(angle <= 0.000001) angle = 0.0000011;
                      bin1 = int((angle-0.000001)*inv_binwidth*MULT_1D_BINS);
                      bin2 = int((angle-0.000001)*inv_binwidth);
                      bin3 = int((angle-0.000001)*SAFE_1D_BINS/360);
                     
                      *(counts  + mybootstrap*nchi*nbins*MULT_1D_BINS + mychi*nbins*MULT_1D_BINS + bin1) += 1;
                      *(counts2 + mybootstrap*nchi*nbins              + mychi*nbins              + bin2) += 1;
                      *(counts3 + mybootstrap*nchi*SAFE_1D_BINS       + mychi*SAFE_1D_BINS       + bin3) += 1;
                      }
                    *(numangles_bootstrap + mybootstrap) += mynumangles;
                    }
                
              }
              """

           #for this next one, counts now will be a properly normalized pdf
           #bin1 = int((angle-0.0000001 + 180)*inv_binwidth*MULT_1D_BINS);
      code_adaptive_doublebins = """
              // weave8b
              int mynumangles, mynumangles_sum, bin1, simnum, found;
              double angle, bin_bottom, bin_top, binwidth;
              for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                mynumangles = 0;
                *(numangles_bootstrap + mybootstrap) = 0 ;   
                for(int mysim=0; mysim < bootstrap_choose; mysim++) {
                    simnum = *(which_runs + mybootstrap*bootstrap_choose + mysim);
                    mynumangles = *(numangles + simnum);   
                    mynumangles_sum = *(numangles_bootstrap + mybootstrap);

                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(angles + mychi*num_sims*max_angles  + simnum*max_angles + anglenum);
                      //angle = 360.0 * (anglenum / mynumangles); 
                      //if(mybootstrap==0 && mysim == 0) printf("angle:%f\\n",angle);
                      if(angle <= 0) angle = 0;
                      if(angle > 360) angle = 360;
                      found = 0;
                      bin_bottom = *(ent_hist_left_breaks + 
                                     + mychi*(nbins*MULT_1D_BINS + 1) + 0);
                      for(int bin1=0; bin1 < nbins*MULT_1D_BINS; bin1++) {
                          bin_top = *(ent_hist_left_breaks + 
                                     + mychi*(nbins*MULT_1D_BINS + 1) + (bin1 + 1));
                          if(found == 0 && angle >= bin_bottom && angle < bin_top) {
                              *(counts_adaptive + mybootstrap*nchi*nbins*MULT_1D_BINS
                                     + mychi*nbins*MULT_1D_BINS + bin1) += 1;
                              found = 1;
                              //if(mybootstrap==0 && mysim == 0 ) printf("bin bot:%f bin top:%f\\n", bin_bottom, bin_top);
                          
                          }
                          bin_bottom = bin_top;
                          if(found == 1) break;
                       }
                              
                    }
                    
                    *(numangles_bootstrap + mybootstrap) += mynumangles;
                }
                for(bin1 = 0; bin1 < nbins*MULT_1D_BINS; bin1++) {
                    binwidth = *(ent_hist_binwidths +  mybootstrap*nchi*nbins*MULT_1D_BINS
                                  + mychi*nbins*MULT_1D_BINS + bin1);
                    *(counts_adaptive + mybootstrap*nchi*nbins*MULT_1D_BINS + mychi*nbins*MULT_1D_BINS + bin1)
                                    /= (*(numangles_bootstrap + mybootstrap) * binwidth );
                    //printf("%f , %f \\n",*(counts_adaptive + mybootstrap*nchi*nbins*MULT_1D_BINS + mychi*nbins*MULT_1D_BINS + bin1),binwidth);
                }
              }
              """

                              
      code_nearest_neighbor_1D = """
          // weave8c2 
          #include <math.h>
              int mynumangles, mynumangles_sum, mynumangles_thisboot, bin1, simnum, found;
              double angle, neighbor_left,neighbor_right,nearest_neighbor;
              double dleft,dright,dmin, logd, dleft_safe, weight;
              double log_nn_dists_sum = 0.0;
              double Lk_minus1 = 0.0;
              for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                mynumangles_thisboot = *(numangles_bootstrap + mybootstrap);
                for(int k_nn = 1; k_nn <= MAX_NEAREST_NEIGHBORS; k_nn++)
                {
                 log_nn_dists_sum = 0;
                 
                 double neighborlist[2*k_nn + 1];
                 double temp_neighbor_dist;
                 double dist_to_neighbors[2*k_nn + 1];
                 
                 for (int anglenum=0; anglenum < mynumangles_thisboot; anglenum++) {

                      angle = *(boot_sorted_angles + mychi*bootstrap_sets*bootstrap_choose*max_angles + mybootstrap*bootstrap_choose*max_angles + anglenum);
                      
                      // make a list of candidate nearest neighbors
                      
                      for(int myk = -k_nn; myk <= k_nn; myk++)
                      {
                       if(myk < 0 && anglenum < abs(myk)) neighborlist[k_nn +myk] = -720;
                       else if(myk > 0 && anglenum > mynumangles_thisboot - myk - 1) neighborlist[k_nn +myk] = 720;
                       else if(myk == 0) neighborlist[k_nn + myk] = 9999;
                       else neighborlist[k_nn + myk] = *(boot_sorted_angles + mychi*bootstrap_sets*bootstrap_choose*max_angles + mybootstrap*bootstrap_choose*max_angles + anglenum + myk);
                       if (myk < 0 )      dist_to_neighbors[k_nn + myk] = angle - neighborlist[k_nn + myk];
                       else if (myk == 0) dist_to_neighbors[k_nn + myk] = 999; // don't want vanishing nn dists
                       else if (myk > 0 ) dist_to_neighbors[k_nn + myk] = neighborlist[k_nn + myk] - angle;
                      }

                      
                      
                      // find the distance to the "k_nn"th nearest neighbor

                      // first, sort using bubble sort
                      // ANN nearest neighbor package uses k-d tree
                      // sort is easier than a k-d tree because it's only 1D :)
                      // http://www.cs.princeton.edu/~ah/alg_anim/gawain-4.0/BubbleSort.html
                      
                      for (int i=0; i< (2*k_nn +1) -1; i++) {
                        for (int j=0; j<(2*k_nn +1) -1 -i; j++)
                           if (dist_to_neighbors[j+1] < dist_to_neighbors[j]) { 
                           temp_neighbor_dist = dist_to_neighbors[j];       
                           dist_to_neighbors[j] = dist_to_neighbors[j+1];
                           dist_to_neighbors[j+1] = temp_neighbor_dist;
                        }
                      }

                     // for(int myk = -k_nn; myk <= k_nn; myk++)
                     // {
                     //    if(anglenum < 10) printf("k:%i dist to neighbor k:%f \\n ",
                     //       k_nn + myk, dist_to_neighbors[k_nn + myk]);
                     // }      

                      dmin = dist_to_neighbors[k_nn - 1 ]; //find kth nearest neighbor
                      if(dmin >= 0.00099) logd = log(dmin);
                      else logd = log(0.00099);
                      //logd = 1000 * log(360.0 / (5005000.0 / ( int((k_nn - 1)/2) + 1)));  //for uniform dist overwrite test
                      log_nn_dists_sum += logd ;

                      //printf("angle:%f nearest:%f next-nearest:%f boot:%i logd:%f, log_nn_dists_sum:%f \\n ",
                      //     angle, dist_to_neighbors[k_nn - 1],dist_to_neighbors[k_nn],mybootstrap, logd, log_nn_dists_sum) ;
                      
                      
                    } // end for(angle ... )
                
                // Now Calculate Entropy From NN Distances For This Bootstrap Sample
                // S = s/n * sum_log_nn_dists + ln (n * vol(s-dim sphere)) -L(k-1)=0 + EULER_GAMMA
                 //L(k-1) = sum(1/i,i=1 .. k)
                 Lk_minus1 = 0.0;
                 if( k_nn > 1) {
                    for(int iter = 1; iter <= k_nn - 1; iter++) {
                           Lk_minus1 += 1.0 / iter; //formula for L(k - 1)
                    }
                 }
                 mynumangles_sum =mynumangles_thisboot;
                 //mynumangles_sum *= 1000; // for uniform dist. overwrite test
                 //printf("log_nn_dists:%f, mynumangles_sum:%i, mybootstrap:%i, mychi:%i, k_nn=%i, Lk_minus1=%f\\n",(float)log_nn_dists_sum, (int)mynumangles_sum, mybootstrap, mychi, k_nn,Lk_minus1);
                  log_nn_dists_sum /= mynumangles_sum;
                  log_nn_dists_sum += log(3.141592653589793238462643383279502884197 / 180.0);
                  log_nn_dists_sum += log( 2 * mynumangles_sum);
                  log_nn_dists_sum += -Lk_minus1 + 0.57721566490153286060651209; //asymptotic bias
                 *(ent_from_sum_log_nn_dists + mybootstrap*nchi*MAX_NEAREST_NEIGHBORS + mychi*MAX_NEAREST_NEIGHBORS + k_nn - 1) = log_nn_dists_sum;
                  log_nn_dists_sum  = 0;
                }
                //*(numangles_bootstrap + mybootstrap) = mynumangles_sum;
              }
              """
      
      
      
      code_nonadaptive_singlebins = """
              // weave9
              int mynumangles, mynumangles_sum, bin1, simnum;
              double angle;
              for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                mynumangles = 0;
                *(numangles_bootstrap + mybootstrap) = 0 ;   
                for (int mysim=0; mysim < bootstrap_choose; mysim++) {
                    simnum = *(which_runs + mybootstrap*bootstrap_choose + mysim);
                    mynumangles = *(numangles + simnum);
                    mynumangles_sum = *(numangles_bootstrap + mybootstrap);
                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(angles + mychi*num_sims*max_angles  + simnum*max_angles + anglenum);
                      if(angle > 360) angle = angle - 360;
                      if(angle <= 0.000001) angle = 0.0000011;
                      bin1 = int((angle-0.000001)*inv_binwidth);
                      if(bin1 < 0)
                      {
                         printf("WARNING: bin less than zero");
                         bin1 = 0;
                      }
                      //printf("bootstrap: %4i boot_sim: %4i sim_num: %4i angle: %3.3f, bin: %4i \\n", mybootstrap, mysim, simnum, angle, bin1);  
                      *(chi_pop_hist + mybootstrap*nchi*nbins + mychi*nbins + bin1) += 1;
                      *(chi_counts + mybootstrap*nchi*nbins + mychi*nbins + bin1) += 1;
                      *(bins + mychi*(permutations + 1)*bootstrap_sets*bootstrap_choose*max_num_angles
                                 + mybootstrap*bootstrap_choose*max_num_angles + mynumangles_sum + anglenum) = bin1;
                      //if(mybootstrap==0 && mysim == 0) printf("angle:%3.3f %3i\\n",angle, bin1);
                      }
                    *(numangles_bootstrap + mybootstrap) += mynumangles;
                  }
                for(bin1 = 0; bin1 < nbins; bin1++) {
                    *(chi_pop_hist + mybootstrap*nchi*nbins + mychi*nbins + bin1) /=  *(numangles_bootstrap + mybootstrap);
                }
              
              }
              """
      #printf("data slot %d bin  %d \\n",temp,bin1);           
      # printf("simnum %d mynumangles %d mynumangles_sum %d\\n",simnum,mynumangles,mynumangles_sum);
      code_adaptive_singlebins = """
              // weave10
              // 
              int mynumangles, mynumangles_sum, bin1, simnum, temp;
              int angle;
              float inv_binwidth_adaptive;
              for(int mybootstrap=0; mybootstrap < bootstrap_sets; mybootstrap++) {
                mynumangles = 0;
                *(numangles_bootstrap + mybootstrap) = 0 ;
                inv_binwidth_adaptive = *(inv_binwidth_adaptive_bootstraps + mybootstrap);
                for (int mysim=0; mysim < bootstrap_choose; mysim++) {
                    simnum = *(which_runs + mybootstrap*bootstrap_choose + mysim);                
                    mynumangles = *(numangles + simnum);
                    mynumangles_sum = *(numangles_bootstrap + mybootstrap);
                    
                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(boot_ranked_angles + mychi*bootstrap_sets*bootstrap_choose*max_angles + mybootstrap*bootstrap_choose*max_angles +mynumangles_sum + anglenum);
                      //printf("chi:%i sim:%i rank:%i\\n",mychi, simnum, angle); 
                      bin1 = int(double(angle)*inv_binwidth_adaptive);
                      if (bin1 < 0) bin1 = 0;
                      else if (bin1 >= nbins) bin1 = nbins - 1;
                      *(chi_pop_hist + mybootstrap*nchi*nbins + mychi*nbins + bin1) += 1;
                      *(chi_counts  + mybootstrap*nchi*nbins + mychi*nbins + bin1) += 1;
                      temp = mychi*(permutations + 1)*bootstrap_sets*bootstrap_choose*max_num_angles
                                 + mybootstrap*bootstrap_choose*max_num_angles + mynumangles_sum + anglenum;

                      *(bins + mychi*(permutations + 1)*bootstrap_sets*bootstrap_choose*max_num_angles
                                 + mybootstrap*bootstrap_choose*max_num_angles + mynumangles_sum + anglenum) = bin1;
                      
                      }
                    *(numangles_bootstrap + mybootstrap) += mynumangles;
                  }
                for(bin1 = 0; bin1 < nbins; bin1++) {
                    *(chi_pop_hist + mybootstrap*nchi*nbins + mychi*nbins + bin1) /=  *(numangles_bootstrap + mybootstrap);
                }
              
              }
              """

      code_nonadaptive_singlebins_sequential = """
              // weave11
              int mynumangles, mynumangles_sum, bin1;
              int nbins_max = nbins_cor*num_histogram_sizes_to_try; // maximum bin size for varying bin size chi counts
              double angle;
                for (int simnum=0; simnum < num_sims; simnum++) {
                    mynumangles = *(numangles + simnum);
                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(angles + mychi*num_sims*max_angles  + simnum*max_angles + anglenum);
                      if(angle > 360) angle = angle - 360;
                      if(angle <= 0.000001) angle = 0.0000011;
                      bin1 = int((angle-0.000001)*(inv_binwidth_cor)); //chi counts for half bin size
                      
                      *(chi_pop_hist_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) += 1;
                      *(chi_counts_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) += 1;
                      *(simbins + mychi*(permutations + 1)*num_sims*max_num_angles
                                 + simnum*max_num_angles + anglenum) = bin1;
                      
                      for (int binmult = 0; binmult < num_histogram_sizes_to_try; binmult++)
                          {
                            bin1 = int((angle-0.000001)*(inv_binwidth_cor*(binmult+1))); //chi counts for half bin size times a multiplication factor
                            *(chi_counts_sequential_varying_bin_size + binmult*num_sims*nchi*nbins_max + simnum*nchi*nbins_max + mychi*nbins_max + bin1) += 1;
                          }
                             
                      
                      }
                    for(bin1 = 0; bin1 < nbins_cor; bin1++) {
                      *(chi_pop_hist_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) /=  mynumangles;
                }
              
              }
              """


      code_adaptive_singlebins_sequential = """
              // weave12b
              int mynumangles, mynumangles_sum, bin1;
              int angle;
                for (int simnum=0; simnum < num_sims; simnum++) {
                    mynumangles = *(numangles + simnum);
                    float inv_binwidth_adaptive_this_sim = *(inv_binwidth_adaptive_sequential + simnum);
                    for (int anglenum=0; anglenum< mynumangles; anglenum++) {
                      angle = *(rank_order_angles_sequential + mychi*num_sims*max_angles  + simnum*max_angles + anglenum);
                      bin1 = int(double(angle)*inv_binwidth_adaptive_this_sim);
                      if (bin1 < 0) bin1 = 0;
                      else if (bin1 >= nbins_cor) bin1 = nbins_cor - 1;
                      *(chi_pop_hist_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) += 1;
                      *(chi_counts_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) += 1;
                      *(simbins + mychi*(permutations + 1)*num_sims*max_num_angles
                                 + simnum*max_num_angles + anglenum) = bin1;
                      
                      }
                   for(bin1 = 0; bin1 < nbins_cor; bin1++) {
                      *(chi_pop_hist_sequential + simnum*nchi*nbins_cor + mychi*nbins_cor + bin1) /=  mynumangles;
                }
              
              }
              """
                                                                              
#
#
#
#
#              weave.inline(code, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins1', 'bins2', 'pop_matrix','permutations','bootstrap_sets'],
#                 #type_converters = converters.blitz,
#                 compiler = mycompiler,runtime_library_dirs="/usr/lib64/", library_dirs="/usr/lib64/", libraries="stdc++")
#
#
#
#
#              weave.inline(code, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins1', 'bins2', 'pop_matrix','permutations','bootstrap_sets'],
#                 #type_converters = converters.blitz,
#                 compiler = mycompiler,runtime_library_dirs="/usr/lib64/", library_dirs="/usr/lib64/", libraries="stdc++")
#
#

      angles = self.angles
      rank_order_angles = self.rank_order_angles
      boot_ranked_angles = self.boot_ranked_angles
      rank_order_angles_sequential = self.rank_order_angles_sequential
      nchi = self.nchi
                                                                              
      counts = self.counts
      counts2 = self.counts2
      counts3 = self.counts3
      chi_counts = self.chi_counts
      numangles = self.numangles
      numangles_bootstrap = self.numangles_bootstrap # this will be overwritten in the weaves,
                                                     #but overwritten correctly before each weave is done
      chi_pop_hist = self.chi_pop_hist
                                                                              
                                                                            
      chi_pop_hist_sequential = self.chi_pop_hist_sequential
      chi_counts_sequential = self.chi_counts_sequential
      chi_counts_sequential_varying_bin_size = self.chi_counts_sequential_varying_bin_size
      bins = self.bins
      simbins = self.simbins
      #ent_hist_left_breaks = self.ent_hist_left_breaks
      #ent_hist_binwidths   = self.ent_hist_binwidths
      boot_sorted_angles = self.boot_sorted_angles
      ent_from_sum_log_nn_dists = self.ent_from_sum_log_nn_dists
      
      for myscramble in range(self.permutations + 1):
          for mychi in range(self.nchi):
             if (myscramble == 1): #resize
                 for bootstrap in range(bootstrap_sets):
                     self.bins[mychi, :, bootstrap, :] = self.bins[mychi, 0, bootstrap, :] #replicate data
                 for sequential_sim_num in range(num_sims):
                      self.simbins[mychi, :, sequential_sim_num, :] = self.simbins[mychi, 0, sequential_sim_num, :]
             if(myscramble > 0):
                 for bootstrap in range(bootstrap_sets): 
                      random.shuffle(self.bins[mychi, myscramble, bootstrap, :])
                 for sequential_sim_num in range(num_sims):
                     random.shuffle(self.simbins[mychi, myscramble, sequential_sim_num, :])
             else:             
                 #for bootstrap in range(bootstrap_sets):
                 #    self.numangles_bootstrap[bootstrap] = 0
                 #    for sequential_sim_num in self.which_runs[bootstrap]:
                 #        for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                 #            bin_num = binsingle(self.angles[mychi,sequential_sim_num,anglenum],MULT_1D_BINS * inv_binwidth) #no adaptive paritioning here
                 #            self.counts[bootstrap, mychi ,bin_num] +=1 #counts are for 1-D histograms for which we use naive binning
                 #        self.numangles_bootstrap[bootstrap] += self.numangles[sequential_sim_num]
                 #    self.counts[bootstrap, mychi, :] /= self.numangles_bootstrap[bootstrap]  # normalize
                
                 weave.inline(code_nonadaptive_doublebins, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins', 'permutations','bootstrap_sets','bootstrap_choose','angles','numangles','counts','counts2','counts3','which_runs','nchi','inv_binwidth','mychi',"max_angles",'max_num_angles','MULT_1D_BINS','FEWER_COR_BTW_BINS','SAFE_1D_BINS','TWOPI'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
                 #weave.inline(code_adaptive_doublebins, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins', 'permutations','bootstrap_sets','bootstrap_choose','angles','numangles','counts_adaptive','which_runs','nchi','ent_hist_binwidths','ent_hist_left_breaks','mychi',"max_angles",'max_num_angles','MULT_1D_BINS'], compiler = mycompiler,runtime_library_dirs="/usr/lib64/", library_dirs="/usr/lib64/", libraries="stdc++")
                 weave.inline(code_nearest_neighbor_1D, ['MAX_NEAREST_NEIGHBORS','numangles_bootstrap','permutations','bootstrap_sets','bootstrap_choose','numangles','which_runs','nchi','mychi',"max_angles",'ent_from_sum_log_nn_dists','boot_sorted_angles'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
                 
                 #print "ent hist stuff compiled successfully"
                 if(adaptive_partitioning == 0):
                    #for bootstrap in range(bootstrap_sets):
                    #       self.numangles_bootstrap[bootstrap] =0
                    #       for sequential_sim_num in self.which_runs[bootstrap]:
                    #           for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                    #               bin_num = binsingle(self.angles[mychi,sequential_sim_num,anglenum],inv_binwidth) #no adaptive paritioning here
                    #               self.chi_pop_hist[bootstrap, mychi , bin_num] +=1 #counts are for 1-D histograms for which we use naive binning
                    #               self.bins[mychi, 0, bootstrap, self.numangles_bootstrap[bootstrap]  + anglenum] = bin_num #use naive binning for 2-D histograms
                    #               self.chi_counts[bootstrap, mychi, bin_num] += 1 
                    #           self.numangles_bootstrap[bootstrap]  += self.numangles[sequential_sim_num]
                    #       self.chi_pop_hist[bootstrap, mychi, :] /=  self.numangles_bootstrap[bootstrap] # normalize
                    
                    weave.inline(code_nonadaptive_singlebins, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins', 'permutations','bootstrap_sets','bootstrap_choose','angles','numangles','chi_pop_hist','chi_counts','which_runs','nchi','inv_binwidth','mychi',"max_angles",'max_num_angles'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
                    
                    # use bins twice as wide for mutinf between sims
                    #       for sequential_sim_num in range(num_sims):
                    #           for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                    #               bin_num = binsingle(self.angles[mychi,sequential_sim_num,anglenum],inv_binwidth/2.0) #no adaptive paritioning here
                    #               self.chi_pop_hist_sequential[sequential_sim_num, mychi , bin_num] +=1 #counts are for 1-D histograms for which we use naive binning
                    #               self.simbins[mychi, 0, sequential_sim_num, anglenum] = bin_num #use naive binning for 2-D histograms
                    #           self.chi_pop_hist_sequential[sequential_sim_num, mychi, :] /=  self.numangles[sequential_sim_num] # normalize
                    weave.inline(code_nonadaptive_singlebins_sequential, ['num_sims', 'numangles_bootstrap', 'nbins_cor', 'simbins', 'permutations','bootstrap_sets','angles','numangles','chi_pop_hist_sequential','chi_counts_sequential','chi_counts_sequential_varying_bin_size','which_runs','nchi','inv_binwidth_cor','mychi',"max_angles",'max_num_angles','FEWER_COR_BTW_BINS','num_histogram_sizes_to_try'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
                 else:
                    # use bins twice as wide for mutinf between sims
                    # for bootstrap in range(bootstrap_sets):
                    #     self.numangles_bootstrap[bootstrap] = 0
                    #     for sequential_sim_num in self.which_runs[bootstrap]:
                    #         for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                    #             bin_num_adaptive = binsingle_adaptive(self.rank_order_angles[mychi,sequential_sim_num,anglenum],inv_binwidth_adaptive)
                    #             if(bin_num_adaptive < 0):
                    #                 print "warning!!: negative bin number wrapped to bin zero"
                    #                 bin_num_adaptive = 0
                    #             if(bin_num_adaptive >= nbins):
                    #                 print "warning!!: bin number overshot nbins, wrapping to bin nbins-1"
                    #                 bin_num_adaptive = nbins - 1
                    #             self.bins[mychi, 0, bootstrap, self.numangles_bootstrap[bootstrap]  + anglenum] = bin_num_adaptive  # overwrite bin value, this is used for adaptive partitioning for 2-D histograms
                    #             self.chi_pop_hist[bootstrap, mychi , bin_num_adaptive] +=1 #counts for 1-D histograms for 2-D mutual information calculation
                    #         self.numangles_bootstrap[bootstrap]  += self.numangles[sequential_sim_num]
                    #     self.chi_pop_hist[bootstrap, mychi , :] /= self.numangles_bootstrap[bootstrap]  # normalize and overwrite pop_hist value, this is used for adaptive partitioning for 2-D histograms
                    weave.inline(code_adaptive_singlebins, ['num_sims', 'numangles_bootstrap', 'nbins', 'bins', 'permutations','bootstrap_sets','bootstrap_choose','boot_ranked_angles','numangles','chi_pop_hist','chi_counts','which_runs','nchi','inv_binwidth_adaptive_bootstraps','mychi',"max_angles",'max_num_angles','FEWER_COR_BTW_BINS'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
                    #   for sequential_sim_num in range(num_sims):
                    #         for anglenum in range(self.numangles[sequential_sim_num]): # numangles per sim is the same regardless of chi or residue
                    #             bin_num_adaptive = binsingle_adaptive(self.rank_order_angles_sequential[mychi,sequential_sim_num,anglenum],inv_binwidth_adaptive_sequential[sequential_sim_num])
                    #             if(bin_num_adaptive < 0):
                    #                 print "warning!!: negative bin number wrapped to bin zero."
                    #                 bin_num_adaptive = 0
                    #             if(bin_num_adaptive >= nbins_cor):
                    #                 print "warning!!: bin number overshot nbins, wrapping to bin nbins - 1. Rank:"+str(bin_num_adaptive/inv_binwidth_adaptive_sequential[sequential_sim_num])
                    #                 bin_num_adaptive = nbins_cor - 1
                    #             self.chi_pop_hist_sequential[sequential_sim_num, mychi , bin_num_adaptive] +=1 #counts are for 1-D histograms using adapative binning
                    #             self.chi_counts_sequential[sequential_sim_num, mychi , bin_num_adaptive] +=1
                    #             self.simbins[mychi, 0, sequential_sim_num, anglenum] = bin_num_adaptive #use adaptive binning for 2-D histograms
                    #         self.chi_pop_hist_sequential[sequential_sim_num, mychi, :] /=  self.numangles[sequential_sim_num] # normalize
                 #   weave.inline(code_nonadaptive_singlebins_sequential, ['num_sims', 'numangles_bootstrap', 'nbins_cor', 'simbins', 'permutations','bootstrap_sets','angles','numangles','chi_pop_hist_sequential','chi_counts_sequential','which_runs','nchi','inv_binwidth_cor','mychi',"max_angles",'max_num_angles'], compiler = mycompiler,runtime_library_dirs="/usr/lib64/", library_dirs="/usr/lib64/", libraries="stdc++")
                 
                 
                    weave.inline(code_adaptive_singlebins_sequential, ['num_sims', 'numangles_bootstrap', 'nbins_cor', 'simbins', 'permutations','bootstrap_sets','rank_order_angles_sequential','chi_counts_sequential','numangles','chi_pop_hist_sequential','which_runs','nchi','inv_binwidth_adaptive_sequential','mychi',"max_angles",'max_num_angles'], compiler = mycompiler,runtime_library_dirs=["/usr/lib64/"], library_dirs=["/usr/lib64/"], libraries=["stdc++"])
             
      
                #print "chi pop hist sequential1:"
                #print chi_pop_hist_sequential[0,0,:]*numangles[0]
                #print "chi pop hist sequential2:"
                #print chi_pop_hist_sequential[1,0,:]*numangles[1]
      #print "numangles_bootstrap:"+str(self.numangles_bootstrap)+"\n"
         # look out for bin values less than zero
      ## Done with permutations
      if len(self.bins[self.bins<0]) > 0:
          print "ERROR: bin values should not be less than zero: ",
          for i in range(num_sims):
              if(VERBOSE >=2):
                  print "Angles\n"
                  print utils.arr1str1( self.angles[mychi, i, :numangles[i]]),
                  print "Bins\n"
                  for duboot in range(bootstrap_sets):
                      print utils.arr1str1( self.bins[mychi, 0, duboot, :])
                  
          sys.exit(1)

      #if(VERBOSE >= 2):
      #    print self.bins[0:2,:,0,:]
          #print counts_adaptive[0,:,:]
          #print "Counts Nonadaptive:\n"
          #print counts

      #calculate entropy for various binwidths, will take the max later
      calc_entropy(self.counts, self.nchi, numangles_bootstrap, calc_variance=calc_variance,entropy=self.entropy,var_ent=self.var_ent,symmetry=self.symmetry)
      calc_entropy(self.counts2, self.nchi, numangles_bootstrap, calc_variance=calc_variance,entropy=self.entropy2,var_ent=self.var_ent,symmetry=self.symmetry)
      calc_entropy(self.counts3, self.nchi, numangles_bootstrap, calc_variance=calc_variance,entropy=self.entropy3,var_ent=self.var_ent,symmetry=self.symmetry)
      #print (sum(counts_adaptive * self.ent_hist_binwidths, axis=2) * numangles_bootstrap)[0,:]
      #calc_entropy_adaptive(counts_adaptive, num_sims, self.nchi, bootstrap_choose, calc_variance=calc_variance,entropy=self.entropy,var_ent=self.var_ent,binwidths=self.ent_hist_binwidths)
      #print "Nearest Neighbor distances k=1 thru 10"
      #print self.ent_from_sum_log_nn_dists[:,:,:]
      #use average of NN estimates for k=7 thru k=10
      #lower values of k are not used as there might be NN dists below that which we can resolve
      #this approach gets around this problem


      #Alright, now compare entropy estimates achieved using different methods
      #For the KNN estimate, I found that the k=4 or k=5 case seemed to be close to the 
      #histogram entropy.
      #So pick the k-value that gives the lowest positive entropy, and then pick the max of the knn estimate
      # or the histogramming estimate.
      # Phi/Psi angles typically are poor for histogramming
      #This should use NN entropy estimate for highly spiked distributions
      #so that they do not have a corrected entropy lower than zero

      print "Entropies for various bin sizes:"
      print self.entropy
      print self.entropy2
      print self.entropy3
      print "K Nearest neighbor entropies:"
      print (ent_from_sum_log_nn_dists)
      for bootstrap in range(bootstrap_sets):
          for mychi in range(nchi):
              thisent1 = self.entropy[bootstrap,mychi]
              thisent2 = self.entropy2[bootstrap,mychi]
              thisent3 = self.entropy3[bootstrap,mychi]
              knn_ent_estimates = self.ent_from_sum_log_nn_dists[bootstrap,mychi,:] 
              #my_knn_ent_estimate = min(knn_ent_estimates[knn_ent_estimates > 0]) 
              #my_knn_ent_estimate = knn_ent_estimates[0]
              ## Use Abel 2009 JCTC eq. 27 weighting for KNN estimates for k=1,2,3
              #my_knn_ent_estimate = sum(WEIGHTS_KNN * knn_ent_estimates)
              my_knn_ent_estimate = knn_ent_estimates
              #if my_knn_ent_estimate > 0:
              #    self.entropy[bootstrap,mychi] = my_knn_ent_estimate
              #else: ## If the KNN algorithm doesn't numerically work well for our data
              self.entropy[bootstrap,mychi] = max((thisent1, thisent2, thisent3))
              if self.entropy[bootstrap,mychi] < 0:
                  print "WARNING: NEGATIVE ENTROPY DETECTED! "


      #self.entropy = entropy
      #self.var_ent += var_ent
      print "Total entropy (before mutinf): %.5f (s2= %.5f )" % (sum(average(self.entropy,axis=0),axis=0), sum(self.var_ent))
      print self.entropy

      

      # derivative of Kullback-Leibler divergence for this residue's chi angles by themselves
      # this has no second-order terms
      # factor of nbins is because sum over index i is retained, while sum over index i just gives a factor of nbins
      numangles_bootstrap_vector = zeros((bootstrap_sets,nbins),float64)
      for bootstrap in range(bootstrap_sets):
          numangles_bootstrap_vector[bootstrap,:]=numangles_bootstrap[bootstrap]
      self.dKLtot_dchis2 = zeros((bootstrap_sets, self.nchi))

      for bootstrap in range(bootstrap_sets):
        for mychi in range(self.nchi):
            nonzero_bins = chi_counts[bootstrap,mychi,:] > 0    
            self.dKLtot_dchis2[bootstrap,mychi] = nbins * sum(((numangles_bootstrap_vector[bootstrap])[nonzero_bins]) * (- 1.0 / ((self.chi_counts[bootstrap,mychi])[nonzero_bins] * 1.0)), axis=0)
      
       
     #free up things we don't need any more, like angles, rank ordered angles, boot angles, etc.
      if(output_timeseries != "yes"):
             del self.angles
      del self.rank_order_angles
      del self.rank_order_angles_sequential
      del self.sorted_angles
      del self.boot_sorted_angles
      del self.boot_ranked_angles
      del self.ent_from_sum_log_nn_dists


#########################################################################################################################################
##### calc_pair_stats: Mutual Information For All Pairs of Torsions ##########################################################
#########################################################################################################################################


# Calculate the mutual information between all pairs of residues.
# The MI between a pair of residues is the sum of the MI between all combinations of res1-chi? and res2-chi?.
# The variance of the MI is the sum of the individual variances.
# Returns mut_info_res_matrix, mut_info_uncert_matrix
def calc_pair_stats(reslist, run_params):
    rp = run_params
    which_runs = rp.which_runs
    #print rp.which_runs
    bootstrap_sets = len(which_runs)
    #initialize the mut info matrix
    check_for_free_mem()
    mut_info_res_matrix = zeros((bootstrap_sets, len(reslist),len(reslist),6,6),float32)
    mut_info_res_matrix_different_sims = zeros((bootstrap_sets, len(reslist),len(reslist),6,6),float32)
    mut_info_uncert_matrix = zeros((bootstrap_sets, len(reslist),len(reslist),6,6),float32)
    dKLtot_dresi_dresj_matrix = zeros((bootstrap_sets, len(reslist),len(reslist)),float32)
    Counts_ij = zeros((rp.nbins,rp.nbins),float64)
    twoD_hist_boot_avg = zeros((len(reslist),len(reslist),6,6,rp.nbins,rp.nbins),float64) #big matrix of 2D populations sum over bootstraps
    #Loop over the residue list
    for res_ind1, myres1 in zip(range(len(reslist)), reslist):
       for res_ind2, myres2 in zip(range(res_ind1, len(reslist)), reslist[res_ind1:]):
        print "\n#### Working on residues %s and %s (%s and %s):" % (myres1.num, myres2.num, myres1.name, myres2.name) , utils.flush()
        max_S = 0.
        if(OFF_DIAG == 1):
          for mychi1 in range(myres1.nchi):
             for mychi2 in range(myres2.nchi):
                 print 
                 print "%s %s , %s %s chi1/chi2: %d/%d" % (myres1.name,myres1.num, myres2.name,myres2.num, mychi1+1,mychi2+1)
                 check_for_free_mem()
                 mutinf_thisdof = var_mi_thisdof = mutinf_thisdof_different_sims = dKLtot_dKL1_dKL2 = 0 #initialize
                 angle_str = ("%s_chi%d-%s_chi%d"%(myres1, mychi1+1, myres2, mychi2+1)).replace(" ","_")
                 print "twoD hist boot avg shape: " + str(twoD_hist_boot_avg.shape ) 
                 if((res_ind1 != res_ind2) or (res_ind1 == res_ind2 and mychi1 > mychi2)):
                     mutinf_thisdof, var_mi_thisdof, mutinf_thisdof_different_sims, dKLtot_dKL1_dKL2, Counts_ij = \
                                 calc_excess_mutinf(myres1.chi_counts[:,mychi1,:],myres2.chi_counts[:,mychi2,:],\
                                                    myres1.bins[mychi1,:,:,:], myres2.bins[mychi2,:,:,:], \
                                                    myres1.chi_counts_sequential[:,mychi1,:],\
                                                    myres2.chi_counts_sequential[:,mychi2,:],\
                                                    myres1.simbins[mychi1,:,:,:], myres2.simbins[mychi2,:,:,:],\
                                                    rp.num_sims, rp.nbins, myres1.numangles_bootstrap,\
                                                    myres1.numangles, rp.sigalpha, rp.permutations,\
                                                    rp.bootstrap_choose, calc_variance=rp.calc_variance,\
                                                    which_runs=rp.which_runs,pair_runs=rp.pair_runs,\
                                                    calc_mutinf_between_sims=rp.calc_mutinf_between_sims, \
                                                    file_prefix=angle_str, plot_2d_histograms=rp.plot_2d_histograms, \
                                                    adaptive_partitioning = rp.adaptive_partitioning)
                     
                 
                 if(res_ind1 == res_ind2 and mychi1 == mychi2):
                     mut_info_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = myres1.entropy[:,mychi1]
                     mut_info_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = myres1.var_ent[:,mychi1]
                     max_S = 0
                     ## still need to calc dKLdiv here
                     dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += 0 #myres1.dKLtot_dchis2[:,mychi1]
                 
                 elif(res_ind1 == res_ind2 and mychi1 > mychi2):
                     mut_info_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = mut_info_res_matrix[:,res_ind1, res_ind2, mychi2, mychi1]
                     mut_info_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = mut_info_uncert_matrix[:,res_ind1, res_ind2, mychi2, mychi1]
                     max_S = 0
                     ## still need to calc dKLdiv here
                     dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += dKLtot_dKL1_dKL2
                     
                 else:
                     S = 0
                     mychi1 = int(mychi1)
                     mychi2 = int(mychi2)
                     blah = myres1.chi_pop_hist[:,mychi1,:]
                     blah = myres2.chi_pop_hist[:,mychi2,:]
                     blah = myres1.bins[mychi1,:,:,:]
                     blah =  myres2.bins[mychi2,:,:,:]
                     blah =  myres1.chi_pop_hist_sequential[:,mychi1,:]
                     dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += dKLtot_dKL1_dKL2
                     dKLtot_dresi_dresj_matrix[:,res_ind2, res_ind1] += dKLtot_dKL1_dKL2 #note res_ind1 neq res_ind2 here
                     mut_info_res_matrix[:,res_ind1 , res_ind2, mychi1, mychi2] = mutinf_thisdof
                     mut_info_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = var_mi_thisdof
                     mut_info_res_matrix_different_sims[:,res_ind1, res_ind2, mychi1, mychi2] = mutinf_thisdof_different_sims
                     twoD_hist_boot_avg[res_ind1, res_ind2, mychi1, mychi2, :, :] = Counts_ij
                     max_S = max([max_S,S])
                     mut_info_res_matrix[:,res_ind2, res_ind1, mychi2, mychi1] = mut_info_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                     #mut_info_uncert_matrix[res_ind1, res_ind2] = mut_info_uncert_matrix[res_ind1, res_ind2]
                     mut_info_uncert_matrix[:,res_ind2, res_ind1, mychi2, mychi1] = mut_info_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                     mut_info_res_matrix_different_sims[:,res_ind2, res_ind1, mychi2, mychi1] = mut_info_res_matrix_different_sims[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                     twoD_hist_boot_avg[res_ind2, res_ind1, mychi2, mychi1, :, :] = swapaxes(twoD_hist_boot_avg[res_ind1, res_ind2, mychi1, mychi2, :, :],0,1) #symmetric matrix
        #print "mutinf=%.3f (uncert=%.3f; max(S)=%.3f" % (average((mut_info_res_matrix[:,res_ind1, res_ind2, : ,:]).flatten()), sum((mut_info_uncert_matrix[0,res_ind1, res_ind2, :, :]).flatten()), max_S),
        #if max_S > 0.26: print "#####",
        print
        

    return mut_info_res_matrix, mut_info_uncert_matrix, mut_info_res_matrix_different_sims, dKLtot_dresi_dresj_matrix, twoD_hist_boot_avg



#########################################################################################################################################
##### Routines for Loading Data #########################################################################################################
#########################################################################################################################################

class ResListEntry:
    name = 'XXX'
    num = 0
    chain = ' '
    def __init__(self,myname,mynum):
        self.name = myname
        self.num = mynum
    def __init__(self,myname,mynum,mychain):
        self.name = myname
        self.num = mynum
        self.chain = mychain

def load_resfile(run_params, load_angles=True, all_angle_info=None):
    rp = run_params
    if rp.num_structs == None: rp.num_structs = 1500000
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
              res_chain = " "
       if load_angles: 
              reslist.append(ResidueChis(res_name,res_num, xvg_resnum, rp.xvg_basedir, rp.num_sims, rp.num_structs, rp.xvgorpdb, rp.binwidth, rp.sigalpha, rp.permutations, rp.phipsi, rp.backbone_only, rp.adaptive_partitioning, rp.which_runs, rp.pair_runs, bootstrap_choose = rp.bootstrap_choose, calc_variance=rp.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp.xvg_chidir, skip=rp.skip,skip_over_steps=rp.skip_over_steps, calc_mutinf_between_sims=rp.calc_mutinf_between_sims,max_num_chis=rp.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp.pdbfile, xtcfile=rp.xtcfile, output_timeseries=rp.output_timeseries))
       else:  reslist.append(ResListEntry(res_name,res_num,res_chain))
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
    parser.add_option("-l", "--load_matrices_numstructs", default = 0, type = "int", help="if you want to load bootstrap matrices from a previous run, give # of structs per sim (yes|no)")
    parser.add_option("--plot_2d_histograms", default = False, action = "store_true", help="makes 2d histograms for all pairs of dihedrals in the first bootstrap")
    parser.add_option("-z", "--zoom_to_step", default = 0, type = "int", help="skips the first n snapshots in xvg files")
    parser.add_option("-m","--max_num_chis", default = 99, type = "int", help="max number of sidechain chi angles per residue or ligand")
    parser.add_option("-f","--pdbfile", default = None, type = "string", help="pdb structure file for additional 3-coord cartesian per residue")
    parser.add_option("-q","--xtcfile", default = None, type = "string", help="gromacs xtc prefix in 'run' subdirectories for additional 3-coord cartesian per residue")
    parser.add_option("-g","--gcc", default = 'intelem', type = "string", help="numpy distutils ccompiler to use. Recommended ones intelem or gcc")
    parser.add_option("-e","--output_timeseries", default = "no", type = "string", help="output corrected dihedral timeseries (requires more memory) yes|no ")
    (options,args)=parser.parse_args()
    mycompiler = options.gcc
    if len(filter(lambda x: x==None, (options.traj_fns, options.xvg_basedir))) != 1:
        parser.error("ERROR exactly one of --traj_fns or --xvg_basedir must be specified")

    print "COMMANDS: ", " ".join(sys.argv)
    

    # Initialize
    resfile_fn = args[0]
    adaptive_partitioning = (options.adaptive == "yes")
    phipsi = 0
    backbone_only = 0
    if options.backbone == "phipsichi":
        phipsi = 2;
        backbone_only =0
    if options.backbone == "phipsi":
        phipsi = 2;
        backbone_only = 1
    if options.backbone == "coarse_phipsi":
        phipsi = -1
        backbone_only = 1
        print "overriding binwidth, using four bins for coarse discretized backbone"
        options.binwidth = 90.0
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
      binwidth=options.binwidth, bins=bins, which_runs=which_runs, xvgorpdb=xvgorpdb, traj_fns=traj_fns, xvg_basedir=options.xvg_basedir, calc_variance=False, xvg_chidir=options.xvg_chidir,bootstrap_choose=options.bootstrap_set_size,pair_runs=pair_runs_array,skip=options.skip,skip_over_steps=options.zoom_to_step,calc_mutinf_between_sims=options.correct_formutinf_between_sims,load_matrices_numstructs=options.load_matrices_numstructs,plot_2d_histograms=options.plot_2d_histograms,max_num_chis=options.max_num_chis,pdbfile=options.pdbfile, xtcfile=options.xtcfile, output_timeseries=options.output_timeseries)

    print run_params


    #====================================================
    #DO ANALYSIS
    #===================================================

    print "Calculating Entropy and Mutual Information"

    
    independent_mutinf_thisdof = zeros((run_params.permutations,1),float32)
    timer = utils.Timer()


    ### load angle data, calculate entropies and mutual informations between residues (and error-propagated variances)
    if run_params.bootstrap_set_size == None:
        if run_params.load_matrices_numstructs == 0: 
               reslist = load_data(run_params)
               print "TIME to load trajectories & calculate intra-residue entropies: ", timer

               prefix = run_params.get_logfile_prefix()
               ##========================================================
               ## Output Timeseries Data in a big matrix
               ##========================================================
               name_num_list = make_name_num_list(reslist)
               if run_params.output_timeseries == "yes":
                      timeseries_chis_matrix = output_timeseries_chis(prefix+"_timeseries",reslist,name_num_list,run_params.num_sims)
                      print "TIME to output timeseries data: ", timer
        
               mut_info_res_matrix, mut_info_uncert_matrix, dKLtot_dresi_dresj_matrix, twoD_hist_boot_avg = calc_pair_stats(reslist, run_params)
               print "TIME to calculate pair stats: ", timer

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

        prefix = run_params.get_logfile_prefix()

        ##========================================================
        ## Output Timeseries Data in a big matrix
        ##========================================================
        name_num_list = make_name_num_list(reslist)
        if run_params.output_timeseries == "yes":
               timeseries_chis_matrix = output_timeseries_chis(prefix+"_timeseries",reslist,name_num_list,run_params.num_sims)
               print "TIME to output timeseries data: ", timer
    
    
        if run_params.load_matrices_numstructs == 0:
            mut_info_res_matrix, mut_info_uncert_matrix, mut_info_res_matrix_different_sims, dKLtot_dresi_dresj_matrix, twoD_hist_boot_avg = calc_pair_stats(reslist, run_params)

        print "TIME to calculate pair stats: ", timer

        # create a master matrix
        #matrix_list += [calc_pair_stats(reslist, run_params)[0]] # add the entropy/mut_inf matrix to the list of matrices
        #bootstraps_mut_inf_res_matrix = zeros(list(matrix_list[0].shape) + [len(matrix_list)], float32)
        #for i in range(len(matrix_list)): bootstraps_mut_inf_res_matrix[:,:,i] = matrix_list[i]

        #mut_info_res_matrix, mut_info_uncert_matrix = bootstraps_mut_inf_res_matrix.mean(axis=2), bootstraps_mut_inf_res_matrix.std(axis=2)
        #prefix = run_params.get_logfile_prefix() + "_sims%s_choose%d" % (",".join(map(str, sorted(runs_superset))), set_size)
        
    ### output results to disk
    
    ##==============================================================
    # setup output or read in previously-calculated mutual information matrices, if given
    ##==============================================================

    name_num_list=[]
    if run_params.load_matrices_numstructs == 0:
        name_num_list = make_name_num_list(reslist)
        #for res in reslist: name_num_list.append(res.name + str(res.num))
    else:
        rownames = []
        colnames = []
        (test_matrix, rownames, colnames) = read_matrix_chis(prefix+"_bootstrap_0_mutinf.txt")
        rownames = colnames
        #print test_matrix
        mut_info_res_matrix = zeros(((len(run_params.which_runs)),test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3]),float32)
        mut_info_res_matrix_different_sims = zeros(((len(run_params.which_runs)),test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3]),float32)
        mut_info_uncert_matrix = zeros(((len(run_params.which_runs)),test_matrix.shape[0],test_matrix.shape[1],test_matrix.shape[2],test_matrix.shape[3]),float32)
        for bootstrap in range(len(run_params.which_runs)):
            (mut_info_res_matrix[bootstrap,:,:,:,:], rownames, colnames) = read_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf.txt")
            (mut_info_res_matrix_different_sims[bootstrap,:,:,:,:], rownames, colnames) = read_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_different_sims.txt")
        name_num_list = rownames

    
    
    ##########################################################################################################   
    ### FINAL STATISTICAL FILTERING USING WILCOXON TEST (NEW) AND OUTPUT MATRICES
    ##########################################################################################################


    #for bootstrap in range(len(run_params.which_runs)):
    if EACH_BOOTSTRAP_MATRIX == 1 and OFF_DIAG == 1 and run_params.load_matrices_numstructs == 0:
        for bootstrap in range(len(run_params.which_runs)):
            output_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf.txt",mut_info_res_matrix[bootstrap],name_num_list,name_num_list)
            output_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_0diag.txt",mut_info_res_matrix[bootstrap],name_num_list,name_num_list, zero_diag=True)
            output_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_different_sims.txt",mut_info_res_matrix_different_sims[bootstrap],name_num_list,name_num_list)
            
           #output_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_uncert.txt",mut_info_uncert_matrix[bootstrap],name_num_list,name_num_list)
    #    if(OUTPUT_DIAG == 1):
    #        output_diag(prefix+"_bootstrap_"+str(bootstrap)+"_entropy_diag_.txt",mut_info_res_matrix[bootstrap],name_num_list)
    #        output_diag(prefix+"_bootstrap_"+str(bootstrap)+"_entropy_diag_uncert.txt",mut_info_uncert_matrix[bootstrap],name_num_list)

    tot_ent = zeros(mut_info_res_matrix.shape[0],float32)
    tot_ent_sig01 = zeros(mut_info_res_matrix.shape[0],float32)
    tot_ent_sig05 = zeros(mut_info_res_matrix.shape[0],float32)
    tot_ent_diag = zeros(mut_info_res_matrix.shape[0],float32)
    mut_info_res_matrix_different_sims_avg = zeros(mut_info_res_matrix.shape[1:],float32)
    mut_info_res_matrix_avg = zeros(mut_info_res_matrix.shape[1:],float32)
    mut_info_res_matrix_sig_01 = zeros(mut_info_res_matrix.shape,float32)
    mut_info_res_matrix_sig_05 = zeros(mut_info_res_matrix.shape,float32)
    mut_info_res_sumoverchis_matrix_sig_01 = zeros(mut_info_res_matrix.shape[0:3],float32)
    mut_info_res_sumoverchis_matrix_sig_05 = zeros(mut_info_res_matrix.shape[0:3],float32)
    mut_info_uncert_matrix_avg = zeros(mut_info_uncert_matrix.shape[1:],float32)
    mut_info_pval_matrix = zeros(mut_info_uncert_matrix.shape[1:],float32)
    mut_info_res_sumoverchis_matrix_sig_avg_01 = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_sumoverchis_matrix_sig_avg_05 = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_sumoverchis_matrix_sig_avg_01_Snorm = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_sumoverchis_matrix_sig_avg_05_Snorm = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_sumoverchis_matrix_avg = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_uncert_sumoverchis_matrix_avg = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_sumoverchis_matrix_different_sims_avg = zeros(mut_info_res_matrix.shape[1:3],float32)
    mut_info_res_maxoverchis_matrix_avg = zeros(mut_info_res_matrix.shape[1:3],float32)
    Kullback_Leibler_local_covar_boots = zeros(mut_info_res_matrix.shape[0:3], float32)
    Kullback_Leibler_local_covar_avg = zeros(mut_info_res_matrix.shape[1:3], float32)

    for bootstrap in range(len(run_params.which_runs)):
        ### invert dKLtot_dresi_dresj_matrix for each bootstrap then average over bootstraps
        ### Don't do inverse here in case matrix is singular -- do it post-processing
        Kullback_Leibler_local_covar_boots[bootstrap,:,:] = dKLtot_dresi_dresj_matrix[bootstrap,:,:]
    Kullback_Leibler_local_covar_avg = average(Kullback_Leibler_local_covar_boots,axis=0)
    
    #print mut_info_res_sumoverchis_matrix_avg.shape
    #print mut_info_res_matrix[:,0,0,0,1]
    testarray = zeros(mut_info_res_matrix.shape[0],float32)
    #print "Applying Wilcoxon test to filter matrix:"
    mutinf_vals = [] # store mutinf values for printing
    for i in range(mut_info_res_matrix.shape[1]):
        for j in range(i, mut_info_res_matrix.shape[2]):
            for k in range(mut_info_res_matrix.shape[3]):
                for m in range(mut_info_res_matrix.shape[4]):
                    if(sum(mut_info_res_matrix[:,i,j,k,m]) > 0):
                        mutinf_boots = mut_info_res_matrix[:,i,j,k,m].copy()
                        mutinf_boots[mutinf_boots < 0] = 0 #use negative and zero values for significance testing but not in the average
                        #zero values of mutinf_boots will include those zeroed out in the permutation test
                        mutinf = mut_info_res_matrix_avg[i,j,k,m] = average(mutinf_boots[mutinf_boots > 0 ],axis=0) 
                        #uncert = mut_info_uncert_matrix_avg[i,j,k,m] = sqrt(cov(mut_info_res_matrix[:,i,j,k,m]) / (run_params.num_sims))
                        #use vfast_cov_for1D_boot as a fast way to calculate a stdev instead of the function std()
                        uncert = pval = None
                        if(mut_info_res_matrix.shape[0] >= 10):
                            uncert = mut_info_uncert_matrix_avg[i,j,k,m] = sqrt((vfast_cov_for1D_boot(reshape(mut_info_res_matrix[:,i,j,k,m],(mut_info_res_matrix.shape[0],1)))[0,0]) / (run_params.num_sims))
                            pval = mut_info_pval_matrix[i,j,k,m] =(stats.wilcoxon(mut_info_res_matrix[:,i,j,k,m]+SMALL))[1] / 2.0 #offset by SMALL to ensure no values of exactly zero will be removed in the test
                            if pval == 0.00 and (i != j): pval = 1.0
                            if pval <= 0.05:
                                mut_info_res_matrix_sig_05[:,i,j,k,m] = mutinf_boots.copy()
                            if pval <= 0.01:
                                mut_info_res_matrix_sig_01[:,i,j,k,m] = mutinf_boots.copy()

                        else:
                            uncert = mut_info_uncert_matrix_avg[i,j,k,m] = sqrt((vfast_cov_for1D_boot(reshape(mut_info_res_matrix[:,i,j,k,m],(mut_info_res_matrix.shape[0],1)))[0,0]) / (run_params.num_sims))
                            pval = 0

                        if pval <= 0.05 and (i != j or k == m):
                            mut_info_res_sumoverchis_matrix_sig_avg_05[i,j] += mutinf
                            mut_info_res_sumoverchis_matrix_sig_05[:,i,j] += mutinf_boots.copy()

                            if(i == j):
                                tot_ent_sig05 += mutinf_boots.copy()
                            else:
                                tot_ent_sig05 -= mutinf_boots.copy()

                        if pval <= 0.05 and i == j and k != m:
                             mut_info_res_sumoverchis_matrix_sig_avg_05[i,j] += -0.5 * mutinf
                             mut_info_res_sumoverchis_matrix_sig_05[:,i,j] += -0.5 * mutinf_boots.copy()
                             tot_ent_sig05  -= 0.5 * mutinf_boots.copy() 


                        if pval <= 0.01 and (i != j or k == m):
                            mut_info_res_sumoverchis_matrix_sig_avg_01[i,j] += mutinf
                            #if(mutinf > 0):
                                #ent1 = 0
                                #ent2 = 0
                                #for count in range(6):
                                #    ent1 += mut_info_res_matrix_avg[i,i,count,count]
                                #    ent2 +=  mut_info_res_matrix_avg[j,j,count,count
                                #mut_info_res_sumoverchis_matrix_sig_avg_01_Snorm[i,j] +=  
                                #mutinf / min(mut_info_res_matrix_avg[i,i,k,k] + mut_info_res_matrix_avg[j,j,m,m])
                            mut_info_res_sumoverchis_matrix_sig_01[:,i,j] += mutinf_boots.copy()

                            if(i == j):
                                tot_ent_sig01 += mutinf_boots.copy()
                            else:
                                tot_ent_sig01 -= mutinf_boots.copy()

                        if pval <= 0.01 and i == j and k != m:
                             mut_info_res_sumoverchis_matrix_sig_avg_01[i,j] += -0.5 * mutinf
                             mut_info_res_sumoverchis_matrix_sig_01[:,i,j] += -0.5 * mutinf_boots.copy()
                             tot_ent_sig01  -= 0.5 * mutinf_boots.copy()

                        if mutinf > 0 and (i != j or k == m):
                            mut_info_res_sumoverchis_matrix_avg[i,j] += mutinf
                            mut_info_res_sumoverchis_matrix_different_sims_avg[i,j] += average(mut_info_res_matrix_different_sims[:,i,j,k,m])
                            mut_info_res_uncert_sumoverchis_matrix_avg[i,j] += uncert ## error bars' sum over torsion pairs

                            if mut_info_res_maxoverchis_matrix_avg[i,j] < mutinf:  mut_info_res_maxoverchis_matrix_avg[i,j] = mutinf

                            if(i == j):
                                tot_ent += mutinf_boots.copy()
                                tot_ent_diag += mutinf_boots.copy()
                            else:
                                tot_ent -= mutinf_boots.copy()
                        if mutinf > 0 and i == j and k != m:
                             mut_info_res_sumoverchis_matrix_avg[i,j] += -0.5 * mutinf
                        #     mut_info_res_sumoverchis_matrix_avg[i,j] += -0.5 * mutinf
                             mut_info_res_sumoverchis_matrix_different_sims_avg[i,j] += -0.5 * average(mut_info_res_matrix_different_sims[:,i,j,k,m])
                             mut_info_res_uncert_sumoverchis_matrix_avg[i,j] += 0.5 * uncert ## error bars' sum over torsion pairs, factor of 0.5 is because of double counting when looping over all pairs i,j in same residue
                             tot_ent -= 0.5 * mutinf_boots.copy()
                             tot_ent_diag -= 0.5 * mutinf_boots.copy()

                        for mat in (mut_info_pval_matrix, mut_info_res_matrix_avg, mut_info_uncert_matrix_avg):
                            if i != j: mat[j,i,m,k] = mat[i,j,k,m] #symmetric matrix

                        if (i==j and k==m): continue
                        elif (i == j and k < m): mutinf_vals.append([mutinf, uncert, pval, "%5s chi%d %5s chi%d (SAME RES)" % (str(name_num_list[i]), k+1, str(name_num_list[j]), m+1)])
                        elif (i != j): mutinf_vals.append([mutinf, uncert, pval, "%5s chi%d %5s chi%d (DIFF RES)" % (str(name_num_list[i]), k+1, str(name_num_list[j]), m+1)])
                    else:
                        mut_info_res_matrix_avg[i,j,k,m] = 0.0
                        mut_info_uncert_matrix_avg[i,j,k,m] = 0.0
                        mut_info_pval_matrix[i,j,k,m] = 0.0
                        #mutinf_vals.append([0, 0, 0, "%5s chi%d %5s chi%d --DEBUG--" % (str(name_num_list[i]), k+1, str(name_num_list[j]), m+1)]) # for debugging
            mut_info_res_sumoverchis_matrix_avg[j,i] = mut_info_res_sumoverchis_matrix_avg[i,j]
            mut_info_res_maxoverchis_matrix_avg[j,i] = mut_info_res_maxoverchis_matrix_avg[i,j]
            mut_info_res_sumoverchis_matrix_different_sims_avg[j,i] = mut_info_res_sumoverchis_matrix_different_sims_avg[i,j]
            mut_info_res_sumoverchis_matrix_sig_avg_01[j,i] = mut_info_res_sumoverchis_matrix_sig_avg_01[i,j]
            mut_info_res_sumoverchis_matrix_sig_avg_05[j,i] = mut_info_res_sumoverchis_matrix_sig_avg_05[i,j]
            mut_info_res_uncert_sumoverchis_matrix_avg[j,i] = mut_info_res_uncert_sumoverchis_matrix_avg[i,j]
            mut_info_res_sumoverchis_matrix_sig_01[:,j,i] = mut_info_res_sumoverchis_matrix_sig_01[:,i,j] 
            mut_info_res_sumoverchis_matrix_sig_05[:,j,i] = mut_info_res_sumoverchis_matrix_sig_05[:,i,j] 

    mut_info_res_sumoverchis_matrix_max_sig01 = zeros(mut_info_res_matrix.shape[1:3],float32)
    for i in range(mut_info_res_matrix.shape[1]):
        for j in range(i, mut_info_res_matrix.shape[2]):
            for chi1 in range(mut_info_res_matrix.shape[3]):
                chis_mi = mut_info_res_matrix_avg[i,j,chi1,:]
                chis_p = mut_info_pval_matrix[i,j,chi1,:]
                chis_mi[chis_p >= 0.01] = 0
                mut_info_res_matrix_avg[i,j,chi1,:] = chis_mi[:] ##OVERWRITING THIS ONE WE DON'T NEED
            mut_info_res_sumoverchis_matrix_max_sig01[i,j] = max( list(mut_info_res_matrix_avg[i,j,:,:].flatten()))

    #for i in range(mut_info_res_sumoverchis_matrix_sig_01.shape[0]): #normalizing mutinf btw res by min. ent of the 2 res.
    #    for j in range(mut_info_res_sumoverchis_matrix_sig_01.shape[0]):
    #        mut_info_res_sumoverchis_matrix_sig_avg_01_Snorm[i,j] = mut_info_res_sumoverchis_matrix_sig_avg_01[i,j] / min(mut_info_res_sumoverchis_matrix_sig_avg_01[i,i], mut_info_res_sumoverchis_matrix_sig_avg_01[j,j])
    #print out total entropies

    output_entropy(prefix+"_entropies_from_bootstraps.txt",tot_ent)
    output_entropy(prefix+"_entropies_from_bootstraps_sig01.txt",tot_ent_sig01)
    output_entropy(prefix+"_entropies_from_bootstraps_sig05.txt",tot_ent_sig05)
    output_entropy(prefix+"_entropies_from_bootstraps_norescorr.txt",tot_ent_diag)

    output_value(prefix+"_entropy.txt",average(tot_ent))
    output_value(prefix+"_entropy_sig01.txt",average(tot_ent_sig01))
    output_value(prefix+"_entropy_sig05.txt",average(tot_ent_sig05))
    output_value(prefix+"_entropy_norescorr.txt",average(tot_ent_diag))
    output_value(prefix+"_entropy_err.txt",sqrt((vfast_cov_for1D_boot(reshape(tot_ent,(tot_ent.shape[0],1)))[0,0]) / (run_params.num_sims)))
    output_value(prefix+"_entropy_sig01_err.txt",sqrt((vfast_cov_for1D_boot(reshape(tot_ent_sig01,(tot_ent_sig01.shape[0],1)))[0,0]) / (run_params.num_sims)))
    output_value(prefix+"_entropy_sig05_err.txt",sqrt((vfast_cov_for1D_boot(reshape(tot_ent_sig05,(tot_ent_sig05.shape[0],1)))[0,0]) / (run_params.num_sims)))
    output_value(prefix+"_entropy_norescorr_err.txt",sqrt((vfast_cov_for1D_boot(reshape(tot_ent_diag,(tot_ent_diag.shape[0],1)))[0,0]) / (run_params.num_sims)))




    # print out mutinf values

    if EACH_BOOTSTRAP_MATRIX == 1 and OFF_DIAG == 1:
        for bootstrap in range(len(run_params.which_runs)):
                    output_matrix_chis(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_sig01.txt",mut_info_res_matrix_sig_01[bootstrap],name_num_list,name_num_list)
                    output_matrix(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_sig05.txt",mut_info_res_matrix_sig_05[bootstrap],name_num_list,name_num_list)
                    output_matrix(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_res_sig01.txt",mut_info_res_sumoverchis_matrix_sig_01[bootstrap],name_num_list,name_num_list)
                    output_matrix(prefix+"_bootstrap_"+str(bootstrap)+"_mutinf_res_sig05.txt",mut_info_res_sumoverchis_matrix_sig_05[bootstrap],name_num_list,name_num_list)
                    #output_matrix(prefix+"_bootstrap_"+str(bootstrap)+"_KLdivpert_boots.txt",Kullback_Leibler_local_covar_boots[bootstrap],name_num_list,name_num_list)
                    





    for mutinf, uncert, pval, name in reversed(sorted(mutinf_vals)): print "BOOTSTRAP DIHEDRAL RESULTS: %s -> mi %.4f   sd/sqrt(n) %.1e   p %.1e" % (name, mutinf, uncert, pval)

    if(OFF_DIAG == 1):
         output_matrix(prefix+"_bootstrap_avg_mutinf_res_sum.txt",            mut_info_res_sumoverchis_matrix_avg ,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_avg_mutinf_res_sum_0diag.txt",      mut_info_res_sumoverchis_matrix_avg ,name_num_list,name_num_list, zero_diag=True)
         output_matrix(prefix+"_bootstrap_avg_mutinf_res_max_0diag.txt",      mut_info_res_maxoverchis_matrix_avg ,name_num_list,name_num_list, zero_diag=True)
         output_matrix(prefix+"_bootstrap_sigavg01_mutinf_res.txt",           mut_info_res_sumoverchis_matrix_sig_avg_01,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_sigavg01_mutinf_res_0diag.txt",     mut_info_res_sumoverchis_matrix_sig_avg_01,name_num_list,name_num_list, zero_diag=True)
         output_matrix(prefix+"_bootstrap_sigavg05_mutinf_res.txt",           mut_info_res_sumoverchis_matrix_sig_avg_05 ,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_sigavg05_mutinf_res_0diag.txt",     mut_info_res_sumoverchis_matrix_sig_avg_05 ,name_num_list,name_num_list, zero_diag=True)
         output_matrix(prefix+"_bootstrap_sigavg01_mutinf_res_norm_0diag.txt", mut_info_res_sumoverchis_matrix_sig_avg_01_Snorm,name_num_list,name_num_list, zero_diag=True)
         output_matrix_chis(prefix+"_bootstrap_avg_mutinf.txt",                    mut_info_res_matrix_avg,name_num_list,name_num_list)
         output_matrix_chis(prefix+"_bootstrap_avg_mutinf_0diag.txt",              mut_info_res_matrix_avg,name_num_list,name_num_list, zero_diag=True)
         output_matrix(prefix+"_bootstrap_sigavg_mutinf_res_uncert.txt",      mut_info_res_uncert_sumoverchis_matrix_avg,name_num_list,name_num_list)
         output_matrix_chis(prefix+"_bootstrap_sigavg_mutinf_pval.txt",       mut_info_pval_matrix,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_sigavg_mutinf_res_max_pval.txt",   amax(amax(mut_info_pval_matrix,axis=-1),axis=-1),name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_avg_mutinf_different_sims.txt",     mut_info_res_sumoverchis_matrix_different_sims_avg,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_sigmax01_mutinf_res_0diag.txt",     mut_info_res_sumoverchis_matrix_max_sig01,name_num_list,name_num_list, zero_diag=True)
         print "twoD hist boot avg shape: " + str(twoD_hist_boot_avg.shape )                                                             
         output_matrix_chis_2dhists(prefix+"_bootstrap_avg_2d_hists.txt",     twoD_hist_boot_avg, name_num_list, name_num_list, nchi=6, nbins = run_params.nbins, zero_diag=True)
         #output_matrix(prefix+"_bootstrap_avg_KLdivpert_res_0diag.txt",     Kullback_Leibler_local_covar_avg,name_num_list,name_num_list, zero_diag=True)
         
    #if(OUTPUT_DIAG == 1):
    #     output_diag(prefix+"_bootstrap_sigavg_entropy_diag_.txt",mut_info_res_matrix_avg,name_num_list)
    #     output_diag(prefix+"_bootstrap_sigavg__entropy_diag_uncert.txt",mut_info_uncert_matrix_avg,name_num_list)

    print "TIME at finish: ", timer
