#!/bin/tcsh

#argument is PDB file
if(-e /software/intel/intel64.csh) source /software/intel/intel64.csh
set mypdb = `echo ${1} | sed -e 's/.pdb//'`
set mut = ${2}
#cat ${mypdb}.pdb | grep CA | sed -e '/NMA/ d' | sed -e '/ACE/ d' | awk '{print NR+1, substr($0,18,3), substr($0,23,4) substr($0,22,1)}' > ${mypdb}_bootstrap10c10.reslist
/raid1/home/mcclendon/scripts/create_reslist.pl ${mypdb}.pdb ${mut}
cp ${mypdb}_reslist.reslist ${mypdb}_bootstrap10c10.reslist
#perl ~cmcclend/scripts/complex_renum.pl temp > ${mypdb}_bootstrap10c10.reslist
if (! -e mutinf_${2}) mkdir mutinf_${2}
cd mutinf_${2}
~mcclendon/usr/bin/python ~mcclendon/trajtools/dihedral_mutent.py -a "yes" -b "phipsichi" -p 0 -n 10 -o 10 -c "yes" ../${mypdb}_bootstrap10c10.reslist -x ../${mut}/ -i 1 -d / -w 30 >& ${mypdb}_mutent_bootstrap10c10_out.txt
cd ..

#END
