#!/bin/tcsh

## Community Analysis Automated Script
## Christopher L. McClendon
## Released under GNU Limited Public License
## 06/27/2014

## Here are the variables you change to tell the script what your filenames are called
## You need a GROMACS-format trajectory file (which you can make from AMBER using VMD) and a structure file and an alignment selection
## Structure file needs to be in the proper view so that when it is opened in pymol, you have the perspective you want without any rotation
## 

set fileprefix=`echo src_ATP_twoMg_yiling_amber_nowat`
set trajectory=`echo ${fileprefix}.trr` #no waters MD trajectory
set struct_nowat=`echo ${fileprefix}2.pdb` #without waters but with Na+ and Cl- present, and Mg+ or other important metal ions                                         
set alignment_selection=`echo r 93-270`  #assumes 1-based index from AMBER, for example the C-lobe of a kinase

echo "done with initial variable assignments"

#set number of threads for MutInf calculation
setenv OMP_NUM_THREADS 5

####### Shouldn't need to change things below this line



### creates an index file based on pdb file to create residue selections for the residues to align, ATP/Mg+, and everything except Na+/Cl- ions
set myindex=`echo ${struct_nowat} | sed -e 's/.pdb/.ndx/'`

cat > create_index.txt <<EOF
${alignment_selection}
r ATP | r Mg+
0 & ! r Na+ & ! r Cl- & ! r WAT & ! r HOH
q
EOF

cat create_index.txt | make_ndx -f ${struct_nowat} -o ${myindex} 



# counts number of entries in index file
set number_of_index_entries=`echo q | make_ndx -f ${struct_nowat} -n ${myindex} |& grep ":" |& tail -n 7 |& head -n 1 | awk '{print $1}'` 
set alignment_index_number=`python -c "print $number_of_index_entries - 2"`
set no_salt_index=`echo ${number_of_index_entries}`

echo "done with creating index file"


### create a pdb file and rotational/translationally fitted trajectory without Na+ or Cl-

set out_pdb=`echo ${fileprefix}3.pdb` #without Na+ or Cl-                                         
sed -e '/Na+/ d' ${struct_nowat} | sed -e '/Cl-/ d' > ${out_pdb}



set out_traj=`echo ${fileprefix}.xtc`


echo ${alignment_index_number} ${no_salt_index}| trjconv -f ${trajectory} -s ${struct_nowat} -n ${myindex} -fit rot+trans -o $out_traj                                                             
echo ${alignment_index_number} ${no_salt_index}| trjconv -f ${trajectory} -s ${struct_nowat} -n ${myindex} -fit rot+trans -o $out_pdb -dump 1                                                      

cat > test_index.txt <<EOF                                                                         
q                                                                                                                      
EOF
                                                                                                                                                          
cat test_index.txt | make_ndx -f $out_pdb -o `echo ${out_pdb} | sed -e 's/.pdb/.ndx/'`

set myindex2=`echo ${out_pdb} | sed -e 's/.pdb/.ndx/'`  #this new index is without Na+ or Cl-
cat > create_index2.txt <<EOF
${alignment_selection}
r ATP | r Mg+
0 & ! r Na+ & ! r Cl- & ! r WAT & ! r HOH
q
EOF
cat create_index2.txt | make_ndx -f ${out_pdb} -o ${myindex2} 

echo "done with rotational/translational superposition and removing Na+ and Cl-"


## identify non-C-alpha atoms to analyze -- called 'hotatoms'
g_select -s ${struct_nowat} -n ${myindex2} -on hotatoms_ca.ndx -select "(name CA) or (resname ALA and name CB) or (resname CYS and name SG) or (resname ASP and name CG) or (resname GLU and name CD) or (resname PHE and name CE1) or (resname GLY and name O) or (resname HIS and name NE2) or (resname HID and name NE2) or (resname HIP and name NE2) or (resname ILE and name CD1) or (resname LYS and name NZ) or (resname LEU and name CD1) or (resname MET and name CE) or (resname ASN and name CG) or (resname PRO and name CG) or (resname GLN and name CD) or (resname ARG and name NH1) or (resname SER and name OG) or (resname THR and name OG1) or (resname VAL and name CG1) or (resname TRP and name CH2) or (resname TYR and name OH) or (resname ATP and name N1) or (resname ATP and name PG) or (resname S2P and name P) or (resname T2P and name P) or (resname CYM and name SG) or (resname HIE and name NE2) or (resname TYR and name P) "

editconf -f ${out_pdb} -n hotatoms_ca.ndx -o hotatoms1g.pdb                                                                                            

sed -e 's/Y2P/TYR/' ${out_pdb} > blah; mv blah ${out_pdb}

echo "done with selecting C-alpha and non-C-alpha atoms to analyze"


### split MD trajectory into six blocks and calculate dihedrals

echo "reading and splitting trajectory into six blocks"

#foreach myfrac (0.125 0.25 0.5 1.0)
foreach myfrac (1.0)
set prefix=`echo dihedrals_${fileprefix}_${myfrac}`
set outprefix=`echo ${fileprefix}_Clobe_align`
set trajectory=$out_traj
#set mylaststep = `g_gmxcheck -f ${trajectory} |& grep Time | tail -n 1 | awk '{print $2}'`
set mylaststep = `gmxcheck -f ${trajectory} |& grep Time | tail -n 1 | awk '{print $2}'`
#set mytime =  `perl -e "print 120 * $mylaststep"`
set mytimestep = `gmxcheck -f ${out_traj} |& grep Time | tail -n 1 | awk '{print $3}'`
set mytime =  `perl -e "print $mytimestep * $mylaststep * $myfrac"`
echo "mylaststep: ", $mylaststep
echo "mytime: ", $mytime
foreach part (1 2 3 4 5 6)
set mybegin = `perl -e "print 1 + ($part - 1) * ($mytime / 6.0);"`
set myend = `perl -e "print $part * ($mytime / 6.0);"`
mkdir ${prefix}
cd ${prefix}
# need to use everything in myindex2 (without Na+ and Cl-) for next step
echo 0 | trjconv -f ../${trajectory} -s ../${out_pdb} -n ../${myindex2} -b $mybegin -e $myend -o ${outprefix}_part${part}.xtc #takes one part of the trajectory
echo 0 | trjconv -s ../${out_pdb} -f ${outprefix}_part${part}.xtc -o ${outprefix}_ca_hotatoms_part${part}.xtc -n ../hotatoms_ca.ndx # creates a trajectory of just the atoms to analyze for community analysis
mkdir run${part}
cd run${part}
rm *.xvg.gz
g_chi -f  ../${outprefix}_part${part}.xtc -s ../../${out_pdb} -all -phi -psi -maxchi 6 -HChi
gzip *.xvg
mv ../${outprefix}_ca_hotatoms_part${part}.xtc ./
cd ..
cd ..

end

end

echo "done splitting trajectory"

## Finished with setup, now on to calculations



## Perform Mutual Information Calculation using MutInf program

perl ~/mutinf/create_reslist_CA_hotatoms.pl ${out_pdb} dihedrals_${fileprefix}_${myfrac} > ${fileprefix}_reslist_out.txt #creates residue list, treating mainchain and sidechain as different resides
set myreslist=`echo ${out_pdb} | sed -e 's/.pdb/_reslist.reslist/'`

## example: of residue list

## 1 ALA 1
## 2 ALA 1S
## 3 MET 2
## 4 MET 2S



#foreach myfrac (0.125 0.25 0.5 1.0)                                                                                                                            
foreach myfrac (1.0)
mkdir mutinf_${fileprefix}_${myfrac}_nonadaptive
cd mutinf_${fileprefix}_${myfrac}_nonadaptive
set prefix=`echo dihedrals_${fileprefix}_${myfrac}`
echo "running MutInf..."
python ~/mutinf/dihedral_mutent.py -x ../${prefix}/ -q ${outprefix}_ca_hotatoms_part -f ../hotatoms1g.pdb -a no -s 0.05 -p 25 -c yes -n 6 -o 6 -w 15 -e no ../${myreslist} -i 1 -d / -g gcc >& mutinf_nonadaptive_out.txt

cd ..

end



## Cluster Mutual Information Matrix using Girvan-Newman Community Analysis



#END



#!/bin/tcsh

#runs in directory with mutinf results

set hotatoms=`echo hotatoms1g.pdb`
#set cutoff=${2} #cutoff distance
#set myfrac=1.0

foreach myfrac (1.0)
#foreach cutoff (16.0 18.0 20.0)
#foreach cutoff (8.0 9.0 10.0)
foreach cutoff (8.0 9.0 10.0)
#foreach cutoff (7.5 8.0 9.0 10.0 12.0 14.0) #(0.125 0.25 0.5 1.0)
   ## Cartesians
   ls mutinf_${fileprefix}_${myfrac}_nonadaptive
   cd mutinf_${fileprefix}_${myfrac}_nonadaptive
   echo $PWD
   set prefix=`ls -lrt *dist_variance_matrix.txt | tail -n 1 | awk '{print $9}' | sed -e 's/_bootstrap_avg_dist_variance_matrix.txt//'`
   echo "prefix"
   echo $prefix
   cd ..
   mkdir  mutinf_${fileprefix}_${cutoff}_${myfrac}_nonadaptive
   cd mutinf_${fileprefix}_${cutoff}_${myfrac}_nonadaptive 

   ## This statement below runs the community analysis  
   echo "running Girvan-Newman community analysis..."
   nohup   python -u ~/mutinf/cluster_list_to_pml.py -f output_clusters_kmeans_${cutoff}.txt -s ../${out_pdb} -m ../mutinf_${fileprefix}_${myfrac}_nonadaptive/${prefix}_bootstrap_avg_dist_variance_matrix.txt -a ../${hotatoms} -t ../mutinf_${fileprefix}_${myfrac}_nonadaptive/${prefix}_bootstrap_avg_mutinf_res_sum_0diag.txt -c ../mutinf_${fileprefix}_${myfrac}_nonadaptive/${prefix}_bootstrap_avg_dist_contacts_cutoff_filter_${cutoff}.txt -o clusters_semirigid1a.pml -u 0.0 -n yes  -x no   >&  nohup_communities_carts.txt &
#
   cd ..

   
#


end

end

#remove GROMACS backup files
rm \#*\#

#END OF COMMUNITY ANALYSIS SCRIPT

#EOF
