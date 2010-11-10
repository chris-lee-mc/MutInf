#!/bin/tcsh

#rm compare_sds_choose3_w1530_p999.txt
#rm compare_ent_choose3_w1530_p999.txt

foreach res (`cat 1M47_diagtest.reslist | sed -f remove_onechis.sed | awk '{print $2 $3}'`) 

set resname = `echo $res | awk '{print substr($1,1,3)}'`
set resnum = `echo $res | awk '{print substr($1,4,3)}'`
echo $resname,$resnum

cat > ${resname}${resnum}.reslist <<EOF
1 $resname $resnum
EOF

rm ${resname}${resnum}_choose3_w1530_p99_out.txt >& duh.null
rm ${resname}${resnum}_choose5_w1530_p999_out.txt >& duh.null

#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 1 1 >> ${resname}${resnum}_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 1 2 >> ${resname}${resnum}_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 1 3 >> ${resname}${resnum}_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 1 4 >> ${resname}${resnum}_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 1 5 >> ${resname}${resnum}_out.txt



#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 2 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 3 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 3 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 3 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 3 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt

python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 3 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 3 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 3 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 3 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 3 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 3 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt

#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 3 4 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 3 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt 
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 3 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt
#python ../../trajtools/dihedral_mutent_nodiag.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 2 3 4 5 >> ${resname}${resnum}_choose3_w1530_p99_out.txt

python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 999 -w 30 ${resname}${resnum}.reslist -n 5 1 2 3 4 5 >> ${resname}${resnum}_choose5_w1530_p999_out.txt




#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 2 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 3 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 1 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 3 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 2 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 3 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 3 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 2 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt

python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 3 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 2 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 3 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 3 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 1 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 3 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 3 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 2 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 3 3 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt

#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 3 4 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 3 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 2 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt 
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 1 3 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt
#python ../../trajtools/dihedral_mutent_verbose.py -x ../../IL2_Flex/1M47_62GLU_SS/ -a "yes" -p 99 -w 30 ${resname}${resnum}.reslist -n 4 2 3 4 5 >> ${resname}${resnum}_choose3_w1530_p999_out.txt




#./postprocess.csh ${resname}${resnum}_choose3_w1530_p99
#./postprocess_ind.csh ${resname}${resnum}_choose3_w1530_p999

end



#END
