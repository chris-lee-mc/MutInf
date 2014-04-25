#Python script to analyze dihedrals. Usage: "python dihedral_analyze.py (list of files to plot). D. Mobley, 5/16/07"
from numpy import *
from matplotlib import *
import sys
import re
from optparse import OptionParser
import os
import gzip
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

from rstyle import rstyle

def ggaxes(fig=None):
    if fig is None: fig = plt.figure()
    ax = fig.add_subplot(111)
    rstyle(ax)
    return ax


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
   if filename.endswith(".xvg.gz"):
       file = gzip.GzipFile(filename, 'r')
   elif filename.endswith(".xvg"):
       file = open(filename,'r');
   #file=open(filename,'r');
   inlines=file.readlines()
   file.close()

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
   dataarray=zeros((inlength,numfields),float64)

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

#===================================================
#READ INPUT ARGUMENTS
#===================================================
parser=OptionParser()
parser.add_option("-p", "--pmf", default="no", type="string", help="whether to use population (no, default) or PMF (yes)")
parser.add_option("-b", "--backbone", default="phipsichi", type="string", help="dihedral angles to use: phipsichi or phispi")
parser.add_option("-d", "--offset", default=0, type=int, help="residue number offset")
parser.add_option("-r", "--reference", default="reference", type="string", help="label for first dataset, i.e. reference")
parser.add_option("-t", "--target", default="target", type="string", help="label for second dataset, i.e. target")
parser.add_option("-s", "--subdir", default="", type="string", help="subdir beneath run1 director for dihedrals")
parser.add_option("-i", "--inputdir", default="run", type="string", help="sum over counts in run directories whose names start with this argument. The default (run) will look for all run directories")
(options,args)=parser.parse_args()


#====================================================
#DO ANALYSIS
#===================================================

#Set binwidth here
binwidth=6. #Degrees
#Calculate number of bins
nbins=int(360./binwidth)

#find out how many directories
directory1 = str(args[0])
directory2 = str(args[1])
residueID = str(args[2])
#print os.listdir(str(directory1)+'/')
dir_names1 = []
dir_names2 = []
#for fn in os.listdir(str(directory1)+'/'):
#    print fn
#    if fn.startswith("run"):
#       dir_names1.append(fn)
#for fn in os.listdir(str(directory2)+'/'):
#    print fn
#    if fn.startswith("run"):
#       dir_names2.append(fn)

dir_names1 = [fn for fn in os.listdir(str(directory1)+'/') if fn.startswith(options.inputdir)]
dir_names2 = [fn for fn in os.listdir(str(directory2)+'/') if fn.startswith(options.inputdir)]
print dir_names1
num_sims=min(len(dir_names1),len(dir_names2))
if (num_sims == len(dir_names1)):
   dir_names = dir_names1
else:
   dir_names = dir_names2

#see what torsions are present

torsion_list = []
if(options.backbone == "phipsichi"):
    torsions_to_check = ["phi","psi","chi1","chi2","chi3","chi4"]
else:
    torsions_to_check = ["phi","psi"]

for torsion in torsions_to_check:
   print "Looking for "+directory1+'/run1/'+options.subdir+torsion+str(residueID)+str(".xvg.gz")
   print os.path.exists(directory1+'/run1/'+options.subdir+torsion+str(residueID)+str(".xvg.gz"))
   print "Looking for "+directory2+'/run1/'+options.subdir+torsion+str(residueID)+str(".xvg.gz")
   print os.path.exists(directory2+'/run1/'+options.subdir+torsion+str(residueID)+str(".xvg.gz"))

   not_appended = True
   for residueID_1 in (residueID, residueID.replace("HIP","HIS"), residueID.replace("HID","HIS"), residueID.replace("HIE","HIS"), residueID.replace("HIS","HID"), residueID.replace("HIS","HIE"), residueID.replace("HIS","HIP"), residueID.replace("T2P","THR"), residueID.replace("S2P","SER"), residueID.replace("THR","T2P"), residueID.replace("SER","S2P")):
      for residueID_2 in (residueID, residueID.replace("HIP","HIS"), residueID.replace("HID","HIS"), residueID.replace("HIE","HIS"), residueID.replace("HIS","HID"), residueID.replace("HIS","HIE"), residueID.replace("HIS","HIP"), residueID.replace("T2P","THR"), residueID.replace("S2P","SER"), residueID.replace("THR","T2P"), residueID.replace("SER","S2P")):
         if((os.path.exists(directory1+'/run1/'+options.subdir+torsion+str(residueID_1)+str(".xvg")) and os.path.exists(directory2+'/run1/'+options.subdir+torsion+str(residueID_2)+str(".xvg"))) or (os.path.exists(directory1+'/run1/'+options.subdir+torsion+str(residueID_1)+str(".xvg.gz")) and os.path.exists(directory2+'/run1/'+options.subdir+torsion+str(residueID_2)+str(".xvg.gz")))):
            residueID_temp = residueID_1
            residueID_targ = residueID_2
            if not_appended: 
                torsion_list.append(torsion)
                not_appended = False

mysearch = re.compile("[0-9]+")
resnum_match=re.search(mysearch, residueID)
residueID.replace("T2P","TPO")
residueID.replace("S2P","SEP")
resnum=str(int(resnum_match.group(0))+int(options.offset))
print "output resnum: "+resnum
residueID_output=re.sub(mysearch,resnum,residueID)


#Figure out stuff about subplots
nplots=len(torsion_list)#+2.0*( torsion_list[0] == "phi" and torsion_list[1] == "psi")
print "torsion list:"
print torsion_list
columns=1 #Fix number of columns at two and add rows
#Compute number of rows
rows=int(ceil(nplots/float(columns)))
print rows


f1=plt.figure(1)
#f1.set_size_inches(4,8.5)
f1.set_size_inches(4,10.5)

ct=0 #To track subplot number
#Loop over all of the input files we want to read

try:
    test = torsion_list[0] #see if we found any torsions
except:
    print "no torsions found for this residue, exiting!"
    sys.exit()

if torsion_list[0] == "phi" and torsion_list[1] == "psi": #then make a ramachandran plot
  ax = plt.subplot(rows, columns, ct+1)

  #Read data
  angles1phi = []
  angles2phi = []
  angles1psi = []
  angles2psi = []
 
  for mydir in dir_names:
     file1phi=directory1+'/'+mydir+'/'+options.subdir+"phi"+str(residueID_temp)+str(".xvg")
     file2phi=directory2+'/'+mydir+'/'+options.subdir+"phi"+str(residueID_targ)+str(".xvg")
     file1psi=directory1+'/'+mydir+'/'+options.subdir+"psi"+str(residueID_temp)+str(".xvg")
     file2psi=directory2+'/'+mydir+'/'+options.subdir+"psi"+str(residueID_targ)+str(".xvg")
     try:
        (data1phi,titlestr)=readxvg(file1phi)
        (data2phi,titlestr2)=readxvg(file2phi)
        (data1psi,titlestr)=readxvg(file1psi)
        (data2psi,titlestr2)=readxvg(file2psi)
     except:
        file1phi=directory1+'/'+mydir+'/'+options.subdir+"phi"+str(residueID_temp)+str(".xvg.gz")
        file2phi=directory2+'/'+mydir+'/'+options.subdir+"phi"+str(residueID_targ)+str(".xvg.gz")
        file1psi=directory1+'/'+mydir+'/'+options.subdir+"psi"+str(residueID_temp)+str(".xvg.gz")
        file2psi=directory2+'/'+mydir+'/'+options.subdir+"psi"+str(residueID_targ)+str(".xvg.gz")
        (data1phi,titlestr)=readxvg(file1phi)
        (data2phi,titlestr2)=readxvg(file2phi)
        (data1psi,titlestr)=readxvg(file1psi)
        (data2psi,titlestr2)=readxvg(file2psi)
     #if angles1phi == None:
     for x in data1phi[:,1]:
         angles1phi.append(x)
     #else:
     #   angles1phi.append(x) for x in list(data1phi[:,1]))

     #if angles2phi == None:
     #   angles2phi = list(data2phi[:,1])
     #else:
     for x in data2phi[:,1]:
        angles2phi.append(x)
     
     #if angles1psi == None:
     #   angles1psi = list(data1psi[:,1])
     #else:
     for x in data1psi[:,1]:
        angles1psi.append(x)

     #if angles2psi == None:
     #   angles2psi = list(data2psi[:,1])
     #else:
     for x in data2psi[:,1]:
        angles2psi.append(x)

  #print "angles1phi"
  #print angles1phi
  temp1=array(angles1phi,float64)
  del angles1phi
  angles1phi = temp1
  temp2=array(angles2phi,float64)
  del angles2phi
  angles2phi = temp2
  temp3=array(angles1psi,float64)
  del angles1psi
  angles1psi = temp3
  temp4=array(angles2psi,float64)
  del angles2psi
  angles2psi = temp4

  #Check and make sure within -180 to 180; I think this should do it 
  angles1phi=(angles1phi+180)%360-180
  angles2phi=(angles2phi+180)%360-180
  angles1psi=(angles1psi+180)%360-180
  angles2psi=(angles2psi+180)%360-180
  
  # ggaxes is a wrapper around rstyle
  #ax = f1.add_subplot(111)
  #ggaxes(plt.figure(figsize=(10,8)))
  #ax.set_xlim(-180,180)
  
  angles1 = zeros((2, len(angles1phi)), float64)
  angles2 = zeros((2, len(angles2phi)), float64)
  
  for i in range(len(angles1phi)):
      angles1[0,i] = angles1phi[i]
      angles1[1,i] = angles1psi[i]
  for i in range(len(angles2phi)):
      angles2[0,i] = angles2phi[i]
      angles2[1,i] = angles2psi[i]
      
  #density1 = gaussian_kde(angles1)
  #density1.covariance_factor = lambda : .01
  #density1._compute_covariance()
  #xs = arange(-180,180,binwidth)
  xs = meshgrid(arange(-180,180,binwidth), arange(-180,180,binwidth))
  #d, m = density1.evaluate(xs)
  density1, xedges, yedges = histogram2d(angles1phi,angles1psi,bins=18, range=((-180,180),(-180,180)), normed=True)
  density1[:,:] *= 1.0
  density1[:,:] += 0.00000001
  #print density1
  if options.pmf == "yes": # -RT ln p  , T=300K
      ys = -0.593 * (300.0/298.0) * log( density1 )# + 0.0000000000001)
  else:
      ys = density1
  
  print "shape of density: "+str(shape(ys))
  print "shape of xedges:  "+str(shape(xedges))
  print "shape of yedges:  "+str(shape(yedges))
  #xlim(-180,180)
  #ylim(-180,180)
  #l1 = ax.contourf(xedges[1:], yedges[1:], ys)
  l1 = ax.imshow(ys, origin='lower', extent=[-180,180,-180,180], aspect='auto') # xedges[1:], yedges[1:], ys)
  #l1 = ax.plot(ys, origin='lower', extent=[-180,180,-180,180]) # xedges[1:], yedges[1:], ys)

  #l1 = ax.plot(xs, ys, antialiased=True, linewidth=2, color="darkgreen") #"#A81450")
  titlestr = str(residueID_output)+" "+" Rama Plot "+options.reference
  #Make plot: A histogram plot
  #n,b,ptc=hist([angles1,angles2],bins=nbins,normed=True,histtype='step')
  #n,b,ptc=hist([density1,density2],bins=nbins,normed=True,histtype='step')
  #Make bar color be black
  #setp(ptc,"facecolor",'k')
  #No y ticks
  plt.setp(plt.gca(), yticklabels=[])

  #Axes labels
  #Do xlabel only on the last row
  #if ct>=nplots-1:
  plt.xlabel('Phi (degrees)')
  plt.ylabel('Psi (degrees)')
  #Minimum/max at 180 degrees
  plt.xlim(-180,180)
  #Labels at 90
  plt.xticks(arange(-180,180.1,90))
  plt.title(titlestr)
  CB = plt.colorbar(l1, orientation='vertical',shrink=1.0)
  #CB.label
  
  ax = plt.subplot(rows, columns, ct+2) 
  #ax.set_xlim(-180,180)
  #density2 = gaussian_kde(angles2)
  #density2.covariance_factor = lambda : .01
  #density2._compute_covariance()
  #xs = arange(-180,180,binwidth)
  xs = meshgrid(arange(-180,180,binwidth), arange(-180,180,binwidth))
  density2, xedges, yedges = histogram2d(angles2phi,angles2psi,bins=18, range=((-180,180),(-180,180)), normed=True)
  density2[:,:] *= 1.0
  density2[:,:] += 0.00000001
  if options.pmf == "yes": # -RT ln p  , T=300K
      ys = -0.593 * (300.0/298.0) * log( density2 )
  else:
      ys = density2
  # ggaxes is a wrapper around rstyle
  #print xs
  #print ys
  l2 = ax.imshow(ys, origin='lower', extent=[-180,180,-180,180],aspect='auto')# xedges[1:], yedges[1:], ys)
  #l2 = ax.plot(ys, origin='lower', extent=[-180,180,-180,180])# xedges[1:], yedges[1:], ys)
  #l2 = ax.contourf(xedges[1:], yedges[1:], ys) 
  
  
  
  #l1 = ax.plot(xs, ys, antialiased=True, linewidth=2, color="#A81450")
  #l1 = ax.fill_between(xs, ys, alpha=.5, zorder=5, antialiased=True, color="#E01B6A")

  

  #http://milkbox.net/smooth_plots_with_matplotlib_pandas/

  #Try and make sense out of the plot title and turn it into something we can use on our plot
  #I am not entirely sure what sort of coding this is using, so I've manually deciphered for chi, psi, omega, phi, etc. This will break if given unexpected coding.
  #Strip quotation marks
  
  #Compile title string
#  titlestr=r'$\rm{1M47  -} \rm{%s }\/%s_%s$' % (res, symbol,subscript)
  titlestr = str(residueID_output)+" "+" Rama Plot "+options.target
  #Make plot: A histogram plot
  #n,b,ptc=hist([angles1,angles2],bins=nbins,normed=True,histtype='step')
  #n,b,ptc=hist([density1,density2],bins=nbins,normed=True,histtype='step')
  #Make bar color be black
  #setp(ptc,"facecolor",'k')
  #No y ticks
  plt.setp(plt.gca(), yticklabels=[])

  #Axes labels
  #Do xlabel only on the last row
  #if ct>=nplots-1:
  plt.xlabel('Phi (degrees)')
  plt.ylabel('Psi (degrees)')
  #Minimum/max at 180 degrees
  plt.xlim(-180,180)
  #Labels at 90
  plt.xticks(arange(-180,180.1,90))
  plt.title(titlestr)
  plt.colorbar(l2, orientation='vertical',shrink=1.0)
  #plt.legend(bbox_to_anchor=(0., 1.00, 1., .100), loc=3,
  #     ncol=2, mode="expand", borderaxespad=0.)
  #Increment subplot counter
  ct+=2

f1.tight_layout()
plotname="rama"+str(residueID_output)+'.eps'
plt.savefig(plotname)
os.system('ps2pdf'+' '+plotname)
os.system('rm '+' '+plotname)

plt.clf()
plt.close()  #close figure

ct = 0
nplots=len(torsion_list)#+2.0*( torsion_list[0] == "phi" and torsion_list[1] == "psi")
print "torsion list:"
print torsion_list
columns=1 #Fix number of columns at two and add rows
#Compute number of rows
rows=int(ceil(nplots/float(columns)))
print rows


f1=plt.figure(1)
#f1.set_size_inches(4,8.5)
f1.set_size_inches(4,10.5)

ct=0 #To track subplot number
#Loop over all of the input files we want to read


    

torsion_index_counter=0
for torsion in torsion_list :
    
  
  ax = plt.subplot(rows, columns, ct+1)

  #Read data
  angles1 = []
  angles2 = []

  
 
  for mydir in dir_names:
     file1=directory1+'/run1/'+options.subdir+torsion+str(residueID_temp)+str(".xvg")
     file2=directory2+'/run1/'+options.subdir+torsion+str(residueID_targ)+str(".xvg")
                      
                      
     try:
        (data1,titlestr)=readxvg(file1)
        (data2,titlestr2)=readxvg(file2)
     except:
        file1=directory1+'/run1/'+options.subdir+torsion+str(residueID_temp)+str(".xvg.gz")
        file2=directory2+'/run1/'+options.subdir+torsion+str(residueID_targ)+str(".xvg.gz")


        (data1,titlestr)=readxvg(file1)
        (data2,titlestr2)=readxvg(file2)
     #if angles1 == None:
     #   angles1 = list(data1[:,1])
     #else:
     for x in data1[:,1]:
        angles1.append(x)
     #if angles2 == None:
     #   angles2 = list(data2[:,1])
     #else:
     for x in data2[:,1]:
        angles2.append(x)
  temp5 = array(angles1, float64)
  del angles1
  angles1 = temp5
  temp6 = array(angles2, float64)
  del angles2
  angles2 = temp6
  #Check and make sure within -180 to 180; I think this should do it 
  angles1=(angles1+180)%360-180
  angles2=(angles2+180)%360-180
  myname = residueID[0:3]
  torsion_index_counter += 1
  print "residue name: "+str(myname)+" torsion index counter :"+str(torsion_index_counter)
  symmetry=False
  if (torsion_index_counter == 4 and (myname == "ASP" or myname == "GLU" or myname == "PHE" or myname == "TYR" \
                      or myname == "NASP" or myname == "NGLU" or myname == "NPHE" or myname == "NTYR" \
                      or myname == "CASP" or myname == "CGLU" or myname == "CPHE" or myname == "CTYR")):
                          #last chi angle is 2-fold symmetric
                          print "correcting for symmetry"
                          angles1[angles1 < -180] = -180
                          angles1[angles1 > 360] = 360
                          angles1[angles1 < 0 ] += 360
                          angles1[angles1 > 180] = angles1[angles1 > 180] - 180
                          angles2[angles2 < -180] = -180
                          angles2[angles2 > 360] = 360
                          angles2[angles2 < 0 ] += 360
                          angles2[angles2 > 180] = angles2[angles2 > 180] - 180
                          symmetry=True
                          
  
  print "angles 1 shape: "+str(shape(angles1))
  angles1len = len(angles1)
  angles2len = len(angles2)
  angles1=resize(angles1,(len(angles1)*3))
  angles2=resize(angles2,(len(angles2)*3))
  print "angles 1 shape: "+str(shape(angles1))
  #wrap around
  for i in range(angles1len):
      angles1[angles1len+i] = angles1[i] - 360
      angles1[angles1len*2+i] = angles1[i] + 360
  for i in range(angles2len):
      angles2[angles2len+i] = angles2[i] - 360
      angles2[angles2len*2+i] = angles2[i] + 360

  
  # ggaxes is a wrapper around rstyle
  #ax = f1.add_subplot(111)
  #ggaxes(plt.figure(figsize=(10,8)))
  #ax.set_xlim(-180,180)
  
  density1 = gaussian_kde(angles1)
  density1.covariance_factor = lambda : .01
  density1._compute_covariance()
  if(symmetry==False):
      xs = arange(-180,180,binwidth)
  else:
      xs = arange(0,180,binwidth)
  if options.pmf == "yes": # -RT ln p  , T=300K
      ys = -0.593 * (300.0/298.0) * log( density1(xs) + 0.0000000000001)
  else:
      ys = density1(xs)
  
  if(symmetry==False):
      plt.xlim(-180,180)
  else:
      plt.xlim(0,180)
  l1 = ax.plot(xs, ys, antialiased=True, linewidth=2, color="darkgreen",label=options.reference) #"#A81450")
  y1hist, edges1 = histogram(angles1, range=[-180,180], bins=int(360/binwidth), density=True)
  y1vals = array(ys, float64)
  print y1hist
  print
  print y1vals

  #ax.set_xlim(-180,180)
  density2 = gaussian_kde(angles2)
  density2.covariance_factor = lambda : .01
  density2._compute_covariance()
  if(symmetry==False):
      xs = arange(-180,180,binwidth)
  else:
      xs = arange(0,180,binwidth)
  if options.pmf == "yes": # -RT ln p  , T=300K
      ys = -0.593 * (300.0/298.0) * log( density2(xs) + 0.0000000000001)
  else:
      ys = density2(xs)
  #ys = density2(xs)
  # ggaxes is a wrapper around rstyle
  #print xs
  #print ys

  l2 = ax.plot(xs, ys, antialiased=True, linewidth=2, color="tan",label=options.target) #"#E01B6A")
  y2hist, edges2 = histogram(angles2, range=[-180,180], bins=int(360/binwidth), density=True)
  y2vals = array(ys, float64)
  print y2hist
  print
  print y2vals
  
  
  flatname1=str(torsion)+residueID_output+'_temp_raw.txt'
  smoothname1=str(torsion)+residueID_output+'_temp_smooth.txt'
  flatname2=str(torsion)+residueID_output+'_targ_raw.txt'
  smoothname2=str(torsion)+residueID_output+'_targ_smooth.txt'

  binfile_res1_flat = open(flatname1,'w')
  binfile_res2_flat = open(flatname2,'w')
  binfile_res1_smooth = open(smoothname1,'w')
  binfile_res2_smooth = open(smoothname2,'w')

  binfile_res1_flat.write(str(nbins)+"\n")
  binfile_res2_flat.write(str(nbins)+"\n")
  binfile_res1_smooth.write(str(nbins)+"\n")
  binfile_res2_smooth.write(str(nbins)+"\n")

  for i in range(nbins):
      binfile_res1_flat.write(str(y1hist[i])+"\n")
      binfile_res2_flat.write(str(y2hist[i])+"\n")
      binfile_res1_smooth.write(str(y1vals[i])+"\n")
      binfile_res2_smooth.write(str(y2vals[i])+"\n")


  binfile_res1_flat.close()
  binfile_res2_flat.close()
  binfile_res1_smooth.close()
  binfile_res2_smooth.close()
    
  
        

        

  #l1 = ax.plot(xs, ys, antialiased=True, linewidth=2, color="#A81450")
  #l1 = ax.fill_between(xs, ys, alpha=.5, zorder=5, antialiased=True, color="#E01B6A")

  

  #http://milkbox.net/smooth_plots_with_matplotlib_pandas/

  #Try and make sense out of the plot title and turn it into something we can use on our plot
  #I am not entirely sure what sort of coding this is using, so I've manually deciphered for chi, psi, omega, phi, etc. This will break if given unexpected coding.
  #Strip quotation marks
  print "title string: "+str(titlestr)
  titlestr=titlestr.replace('"','')
  #Break into the prefix, plus the residue
  (prefix,res)=titlestr.split()
  prefix=prefix.replace('\\f{}',' ') #Remove garbage
  #Figure out symbol
  if prefix[2]=='c': 
	symbol='\chi'
	pname='chi'
  elif prefix[2]=='w': 
        symbol='\omega'
	pname='omega'
  elif prefix[2]=='y': 
	symbol='\psi'
	pname='psi'
  elif prefix[2]=='f': 
	symbol='\phi'
	pname='phi'
  #Get rid of symbol-related stuff
  prefix=prefix[3:]
  #Figure out if there should be a subscript
  try:
    if prefix[1]=='s':
      subscript=prefix[2]
    else: subscript=''
  except: subscript=''

  #Compile title string
#  titlestr=r'$\rm{1M47  -} \rm{%s }\/%s_%s$' % (res, symbol,subscript)
  titlestr = str(residueID_output)+" "+str(torsion)
  #Make plot: A histogram plot
  #n,b,ptc=hist([angles1,angles2],bins=nbins,normed=True,histtype='step')
  #n,b,ptc=hist([density1,density2],bins=nbins,normed=True,histtype='step')
  #Make bar color be black
  #setp(ptc,"facecolor",'k')
  

  #Axes labels
  #Do xlabel only on the last row
  if ct>=nplots-1:
    plt.xlabel('Angle (degrees)')
  if(options.pmf == "no"):
      #No y ticks
      plt.setp(plt.gca(), yticklabels=[])
      plt.ylabel('Probability')
  else:
      yrange = max(ys) - min(ys)
      plt.ylim(min(ys), max(ys))
      plt.yticks(arange(min(ys),max(ys),yrange))
      plt.ylabel('Potential of Mean Force, kcal/mol')
  #Minimum/max at 180 degrees
  if(symmetry==False):
      plt.xlim(-180,180)
      plt.xticks(arange(-180,180.1,90))
  else:
      plt.xlim(0, 180)
      plt.xticks(arange(0,180.1,45))
  #Labels at 90
  
  
  
  
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.,title=titlestr,fontsize=10 )

  #plt.title(titlestr)
  #ax.legend(loc='upper center', shadow=False)
  #plt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0)
  #Increment subplot counter
  ct+=1


#print 'pref',prefix
f1.tight_layout()


left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.75   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
#plotname=pname+res+'.eps'
#plotname=pname+residueID_output+'.eps'
#just call it chi even if the residue only has phi and psi
plotname="chi"+residueID_output+'.eps'

plt.savefig(plotname)
os.system('ps2pdf'+' '+plotname)
os.system('rm '+' '+plotname)
#show()


########################################################################
# GUI
# modified from Ringer
############################3

#from wxtbx import plots
#import wx

# def run_app (results) :
#   app = wx.App(0)
#   frame = RingerFrame(None, -1, "Ringer results")
#   frame.show_results(results)
#   frame.Show()
#   app.MainLoop()

# class RingerFrame (plots.plot_frame) :
#   def create_plot_panel (self) :
#     plot = RingerPlot(self, figure_size=(6,8))
#     plot.canvas.Bind(wx.EVT_CHAR, self.OnChar)
#     return plot

#   def draw_top_panel (self) :
#     self.top_panel = wx.Panel(self, style=wx.SUNKEN_BORDER)
#     panel_szr = wx.BoxSizer(wx.VERTICAL)
#     self.top_panel.SetSizer(panel_szr)
#     szr2 = wx.BoxSizer(wx.HORIZONTAL)
#     panel_szr.Add(szr2)
#     txt1 = wx.StaticText(self.top_panel, -1, "Residue to display:")
#     szr2.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#     self.chooser = wx.Choice(self.top_panel, -1, size=(200,-1))
#     szr2.Add(self.chooser, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#     self.Bind(wx.EVT_CHOICE, self.OnSelect, self.chooser)
#     self.Bind(wx.EVT_CHAR, self.OnChar)
#     self.chooser.Bind(wx.EVT_CHAR, self.OnChar)
#     return self.top_panel

#   def OnSelect (self, event) :
#     selection = event.GetEventObject().GetSelection()
#     self.plot_panel.show_residue(self.results[selection])

#   def show_results (self, results) :
#     self.results = results
#     choices = [ result.format() for result in results ]
#     self.chooser.SetItems(choices)
#     self.chooser.SetSelection(0)
#     self.plot_panel.show_residue(self.results[0])

#   def OnChar (self, event) :
#     key = event.GetKeyCode()
#     if (len(self.results) == 0) : return
#     selection = self.chooser.GetSelection()
#     if (key in [wx.WXK_TAB, wx.WXK_RETURN, wx.WXK_SPACE]) :
#       if (selection < (len(self.results) - 1)) :
#         selection += 1
#       elif (len(self.results) > 0) :
#         selection = 0
#     elif (key in [wx.WXK_DELETE, wx.WXK_BACK]) :
#       if (selection > 0) :
#         selection -= 1
#       else :
#         selection = len(results) - 1
#     self.chooser.SetSelection(selection)
#     self.plot_panel.show_residue(self.results[selection])

# class RingerPlot (plots.plot_container) :
#   def show_residue (self, residue) :
#     if (self.disabled) : return
#     self.figure.clear()
#     subplots = []
#     for i in range(1, residue.n_chi + 1) :
#       chi = residue.get_angle(i)
#       if (chi is None) : continue
#       if (len(subplots) > 0) :
#         p = self.figure.add_subplot(4, 1, i, sharex=subplots[0])
#       else :
#         p = self.figure.add_subplot(4, 1, i)
#         p.set_title(residue.format())
#       p.set_position([0.15, 0.725 - 0.225*(i-1), 0.8, 0.225])
#       x = [ k*chi.sampling for k in range(len(chi.densities)) ]
#       p.plot(x, chi.densities, 'r-', linewidth=1)
#       p.axvline(chi.angle_current, color='b', linewidth=2, linestyle='--')
#       p.axhline(0, color=(0.4,0.4,0.4), linestyle='--', linewidth=1)
#       p.axhspan(0.3,1,facecolor="green",alpha=0.5)
#       p.axhspan(-1,0.3,facecolor="grey",alpha=0.5)
#       p.set_xlim(0,360)
#       ax = p.get_axes()
#       ax.set_ylabel("Rho")
#       ax.set_xlabel("Chi%d" % i)
#       subplots.append(p)
#     for p in subplots[:-1] :
#       for label in p.get_axes().get_xticklabels() :
#         label.set_visible(False)
#     self.canvas.draw()
#     self.canvas.Fit()
#     self.Layout()
#     self.parent.Refresh()

# if (__name__ == "__main__") :
#   run(sys.argv[1:])
