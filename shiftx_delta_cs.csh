#!/bin/tcsh

set pdb1=${1}
set pdb2=${2}
set mydir1=${3}
set mydir2=${4}

set outprefix1=`echo ${pdb1} | sed -e 's/.pdb$//'`
set outprefix2=`echo ${pdb2} | sed -e 's/.pdb$//'`

cd ${mydir1}
cd shiftx2

set cs=`ls -rt *.pdb.cs | tail -n 1`
set reslist1=`echo ${pdb1} | sed -e 's/.pdb$/_reslist.reslist/'`

~/mutinf/create_reslist.pl ../../${pdb1} dihedrals_${mydir1}

egrep "N\,[0-9]+" ${cs} | awk -F "," '{print $4}' > shifts.txt
awk '{print $2 $3}' ../$reslist1 > res.txt
paste res.txt shifts.txt > shifts_res.txt

~/mutinf/replace_res_bfac_mut.pl ../../${pdb1} shifts_res.txt

cp `echo ../../${pdb1} |sed -e 's/.pdb/_resbfac.pdb/'` ../../`echo ${pdb1} | sed -e 's/.pdb/_shiftx_resbfac.pdb/'`

cd ..

cd ..

cd ${mydir2}
cd shiftx2

set cs=`ls -rt *.pdb.cs | tail -n 1`
echo "chemical shifts file "
echo $cs
head -n 10 $cs

set reslist2=`echo ${pdb2} | sed -e 's/.pdb$/_reslist.reslist/'`

~/mutinf/create_reslist.pl ../../${pdb2} dihedrals_${mydir2}

egrep "N\,[0-9]+" ${cs} | awk -F "," '{print $4}' > shifts.txt
awk '{print $2 $3}' ../$reslist2 > res.txt
paste res.txt shifts.txt > shifts_res.txt

~/mutinf/replace_res_bfac_mut.pl ../../${pdb2} shifts_res.txt

cp `echo ../../${pdb2} |sed -e 's/.pdb/_resbfac.pdb/'` ../../`echo ${pdb2} | sed -e 's/.pdb/_shiftx_resbfac.pdb/'`

cd ..
cd ..

perl ~/mutinf/replace_res_diffbfac_mut_cs.pl ${pdb1} ${mydir1}/shiftx2/shifts_res.txt ${mydir2}/shiftx2/shifts_res.txt BLA00 ${outprefix1}_${outprefix2}_delta_cs


#END
