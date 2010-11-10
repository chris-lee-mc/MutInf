#!/bin/tcsh

rm compare_sds_choose3_w1530_p99.txt
rm compare_ent_choose3_w1530_p99.txt


#@ ent_sum_1 = 0 # $ent_sum_1 + `echo $file | head -n 1 | tail -n 1`
#@ ent_sum_2 = 0 # $ent_sum_2 + `echo $file | head -n 2 | tail -n 1`
#@ ent_sum_3 = 0 # $ent_sum_3 +`echo $file | head -n 3 | tail -n 1`
#@ ent_sum_4 = 0 # $ent_sum_4 + `echo $file | head -n 4 | tail -n 1`
#@ ent_sum_5 = 0 # $ent_sum_5 + `echo $file | head -n 5 | tail -n 1`
#@ ent_sum_6 = 0 # $ent_sum_6 + `echo $file | head -n 6 | tail -n 1`
#@ ent_sum_7 = 0 # $ent_sum_7 + `echo $file | head -n 7 | tail -n 1`
#@ ent_sum_8 = 0 # $ent_sum_8 + `echo $file | head -n 8 | tail -n 1`
#@ ent_sum_9 = 0 #$ent_sum_9 + `echo $file | head -n 9 | tail -n 1`
#@ ent_sum_10 = 0 #  $ent_sum_10 + `echo $file | head -n 10 | tail -n 1`
@ ent_sum_11 = 0 # $ent_sum_11 + `echo $file | head -n 11 | tail -n 1`
@ ent_sum_12 = 0 # $ent_sum_12 + `echo $file | head -n 12 | tail -n 1`
@ ent_sum_13 = 0 # $ent_sum_13 + `echo $file | head -n 13 | tail -n 1`
@ ent_sum_14 = 0 # $ent_sum_14 + `echo $file | head -n 14 | tail -n 1`
@ ent_sum_15 = 0 # $ent_sum_15 + `echo $file | head -n 15 | tail -n 1`
@ ent_sum_16 = 0 # $ent_sum_16 + `echo $file | head -n 16 | tail -n 1`
@ ent_sum_17 = 0 # $ent_sum_17 + `echo $file | head -n 17 | tail -n 1`
@ ent_sum_18 = 0 # $ent_sum_18 + `echo $file | head -n 18 | tail -n 1`
@ ent_sum_19 = 0 # $ent_sum_19 + `echo $file | head -n 19 | tail -n 1`
@ ent_sum_20 = 0 #$ent_sum_20 + `echo $file | head -n 20 | tail -n 1`
#@ ent_sum_21 = 0 # $ent_sum_21 + `echo $file | head -n 21 | tail -n 1`
#@ ent_sum_22 = 0 # $ent_sum_22 + `echo $file | head -n 22 | tail -n 1`
#@ ent_sum_23 = 0 # $ent_sum_23 + `echo $file | head -n 23 | tail -n 1`
#@ ent_sum_24 = 0 # $ent_sum_24 + `echo $file | head -n 24 | tail -n 1`
#@ ent_sum_25 = 0 # $ent_sum_25 + `echo $file | head -n 25 | tail -n 1`


foreach res (`cat 1M47_diagtest.reslist | sed -f remove_onechis.sed | awk '{print $2 $3}'`) 

#foreach res (LYP9)
set resname = `echo $res | awk '{print substr($1,1,3)}'`
set resnum = `echo $res | awk '{print substr($1,4,3)}'`
echo $resname,$resnum

#./postprocess.csh ${resname}${resnum}_choose3_w1530_p99
./postprocess_exc.csh  ${resname}${resnum}_choose3_w1530_p99
#./postprocess_ind.csh ${resname}${resnum}_choose5_w1530_p999
./postprocess_entropy.csh ${resname}${resnum}_choose3_w1530_p99
#end



#foreach res (`cat 1M47_diagtest.reslist | sed -f remove_onechis.sed | awk '{print $2 $3}'`) 
#foreach res (LYP9)
set resname = `echo $res | awk '{print substr($1,1,3)}'`
set resnum = `echo $res | awk '{print substr($1,4,3)}'`
echo $resname,$resnum
set file = ${resname}${resnum}_choose3_w1530_p99_entropy.txt
set mean_ent = `echo mean\(read.table\(\"${file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set sd_ent = `echo sd\(read.table\(\"${file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set pvalue_shapiro_ent = `echo shapiro.test\(as.matrix\(read.table\(\"${file}\"\)\)\) | R --no-save | grep value | awk '{print $6}'`
echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(shape\=$mean_ent \* $mean_ent \/ \($sd_ent \* $sd_ent\)\, scale\=$sd_ent \* $sd_ent \/ $mean_ent\)\) 
set pvalue_supremum_anderson_darling_ent = `echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(shape\=$mean_ent \* $mean_ent \/ \($sd_ent \* $sd_ent\)\, scale\=$sd_ent \* $sd_ent \/ $mean_ent\)\) | R --no-save | grep value | awk '{print $6}'`
set mean_ent_agg = `cat ${resname}${resnum}_choose5_w1530_p999_out.txt | grep "Entropy" | awk '{print $2}' | tail -n 1`
set sd_ent_agg = `cat ${resname}${resnum}_choose5_w1530_p999_out.txt | grep "Entropy" | awk '{print $4}' | tail -n 1`
echo $file, $mean_ent, $sd_ent,$mean_ent_agg,$sd_ent_agg,$pvalue_shapiro_ent,$pvalue_supremum_anderson_darling_ent 
echo $file, $mean_ent, $sd_ent,$mean_ent_agg,$sd_ent_agg,$pvalue_shapiro_ent,$pvalue_supremum_anderson_darling_ent  >> compare_ent_choose3_w1530_p99.txt
#set ent_add_1 =  `echo $file | head -n 1 | tail -n 1`
#set ent_add_2 =  `echo $file | head -n 2 | tail -n 1`
#set ent_add_3 = `echo $file | head -n 3 | tail -n 1`
#set ent_add_4 =  `echo $file | head -n 4 | tail -n 1`
#set ent_add_5 =  `echo $file | head -n 5 | tail -n 1`
#set ent_add_6 =  `echo $file | head -n 6 | tail -n 1`
#set ent_add_7 =  `echo $file | head -n 7 | tail -n 1`
#set ent_add_8 =  `echo $file | head -n 8 | tail -n 1`
#set ent_add_9 =  `echo $file | head -n 9 | tail -n 1`
set ent_add_10 =  `echo $file | head -n 10 | tail -n 1`
set ent_add_11 =  `echo $file | head -n 11 | tail -n 1`
set ent_add_12 =  `echo $file | head -n 12 | tail -n 1`
set ent_add_13 =  `echo $file | head -n 13 | tail -n 1`
set ent_add_14 =  `echo $file | head -n 14 | tail -n 1`
set ent_add_15 =  `echo $file | head -n 15 | tail -n 1`
set ent_add_16 = `echo $file | head -n 16 | tail -n 1`
set ent_add_17 =  `echo $file | head -n 17 | tail -n 1`
set ent_add_18 =  `echo $file | head -n 18 | tail -n 1`
set ent_add_19 =  `echo $file | head -n 19 | tail -n 1`
set ent_add_20 =  `echo $file | head -n 20 | tail -n 1`
#set ent_add_21 =  `echo $file | head -n 21 | tail -n 1`
#set ent_add_22 =  `echo $file | head -n 22 | tail -n 1`
#set ent_add_23 =  `echo $file | head -n 23 | tail -n 1`
#set ent_add_24 =  `echo $file | head -n 24 | tail -n 1`
#set ent_add_25 =  `echo $file | head -n 25 | tail -n 1`

#set ent_sum_1 = `perl -e "print ent_sum_1 + $ent_add_1"`
#set ent_sum_2 = `perl -e "print ent_sum_2 + $ent_add_2"`
#set ent_sum_3 = `perl -e "print ent_sum_3 + $ent_add_3"`
#set ent_sum_4 = `perl -e "print ent_sum_4 + $ent_add_4"`
#set ent_sum_5 = `perl -e "print ent_sum_5 + $ent_add_5"`
#set ent_sum_6 = `perl -e "print ent_sum_6 + $ent_add_6"`
#set ent_sum_7 = `perl -e "print ent_sum_7 + $ent_add_7"`
#set ent_sum_8 = `perl -e "print ent_sum_8 + $ent_add_8"`
#set ent_sum_9 = `perl -e "print ent_sum_9 + $ent_add_9"`
# set ent_sum_10 = `perl -e "print ent_sum_10 + $ent_add_10"`
 set ent_sum_11 = `perl -e "print ent_sum_11 + $ent_add_11"`
 set ent_sum_12 = `perl -e "print ent_sum_12 + $ent_add_12"`
 set ent_sum_13 = `perl -e "print ent_sum_13 + $ent_add_13"`
 set ent_sum_14 = `perl -e "print ent_sum_14 + $ent_add_14"`
 set ent_sum_15 = `perl -e "print ent_sum_15 + $ent_add_15"`
 set ent_sum_16 = `perl -e "print ent_sum_16 + $ent_add_16"`
 set ent_sum_17 = `perl -e "print ent_sum_17 + $ent_add_17"`
 set ent_sum_18 = `perl -e "print ent_sum_18 + $ent_add_18"`
 set ent_sum_19 = `perl -e "print ent_sum_19 + $ent_add_19"`
 set ent_sum_20 = `perl -e "print ent_sum_20 + $ent_add_20"`
# set ent_sum_21 = `perl -e "print ent_sum_21 + $ent_add_21"`
# set ent_sum_22 = `perl -e "print ent_sum_22 + $ent_add_22"`
# set ent_sum_23 = `perl -e "print ent_sum_23 + $ent_add_23"`
# set ent_sum_24 = `perl -e "print ent_sum_24 + $ent_add_24"`
# set ent_sum_25 = `perl -e "print ent_sum_25 + $ent_add_25"`

set pdfoutfile = `echo $file | sed -e 's/.txt/_hist.pdf/'`
rm $pdfoutfile
echo pdf\(\"${pdfoutfile}\"\)\;hist\(t\(as.matrix\(read.table\(\"${file}\"\)\)\)\,main\=\"M.I. for Subsamples of ${file}\"\,xlab\=\"Mutual Information \(nats\)\" \)\;dev.off\(\) | R --no-save | tail -n 2 | head -n 1


foreach file (${resname}${resnum}_choose3_w1530_p99_*_*_mut.txt)
#echo sd\(read.table\(\"${file}\"\)\) 
set ind_output_file = `echo $file | sed -e 's/_mut.txt/_ind.txt/' | sed -e 's/choose3/choose5/' | sed -e 's/p99/p999/'`
set exc_output_file = `echo $file | sed -e 's/_mut.txt/_exc.txt/'`
echo $ind_output_file
set sd_bootstrap = `echo sd\(read.table\(\"${file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set mean_bootstrap = `echo mean\(read.table\(\"${file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set pvalue_shapiro = `echo shapiro.test\(as.matrix\(read.table\(\"${file}\"\)\)\) | R --no-save | grep value | awk '{print $6}'`
echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(scale\=$mean_bootstrap \* $mean_bootstrap \/ \($sd_bootstrap \* $sd_bootstrap\)\, shape\=$sd_bootstrap \* $sd_bootstrap \/ $mean_bootstrap\) 
set pvalue_supremum_anderson_darling = `echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(shape\=$mean_bootstrap \* $mean_bootstrap \/ \($sd_bootstrap \* $sd_bootstrap\)\, scale\=$sd_bootstrap \* $sd_bootstrap \/ $mean_bootstrap\)\) | R --no-save | grep value | awk '{print $6}'`
set sd_bootstrap_exc = `echo sd\(read.table\(\"${exc_output_file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set mean_bootstrap_exc = `echo mean\(read.table\(\"${exc_output_file}\"\)\) | R --no-save | tail -n 2 | head -n 1`
set pvalue_shapiro_exc = `echo shapiro.test\(as.matrix\(read.table\(\"${exc_output_file}\"\)\)\) | R --no-save | grep value | awk '{print $6}'`
echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(scale\=$mean_bootstrap \* $mean_bootstrap \/ \($sd_bootstrap \* $sd_bootstrap\)\, shape\=$sd_bootstrap \* $sd_bootstrap \/ $mean_bootstrap\) 
set pvalue_supremum_anderson_darling_exc = `echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${exc_output_file}\"\)\)\,\"dgamma\"\,list\(shape\=$mean_bootstrap_exc \* $mean_bootstrap_exc \/ \($sd_bootstrap_exc \* $sd_bootstrap_exc\)\, scale\=$sd_bootstrap_exc \* $sd_bootstrap_exc \/ $mean_bootstrap_exc\)\) | R --no-save | grep value | awk '{print $6}'`
echo "Bootstrap:",$mean_bootstrap, $sd_bootstrap, $pvalue_shapiro,$pvalue_supremum_anderson_darling
echo "Excess:",$mean_bootstrap_exc, $sd_bootstrap_exc, $pvalue_shapiro_exc,$pvalue_supremum_anderson_darling_exc
set sd_ind = `echo sd\(t\(as.matrix\(read.table\(\"${ind_output_file}\"\)\)\)\) | R --no-save | tail -n 2 | head -n 1`
set mean_ind = `echo mean\(t\(as.matrix\(read.table\(\"${ind_output_file}\"\)\)\)\) | R --no-save | tail -n 2 | head -n 1 | awk '{print $2}'`
set pvalue_shapiro_ind = `echo shapiro.test\(t\(as.matrix\(read.table\(\"${ind_output_file}\"\)\)\)\) | R --no-save | grep value | awk '{print $6}'`
set pvalue_supremum_anderson_darling_ind = `echo library\(\'truncgof\'\)\;adup.test\(as.matrix\(read.table\(\"${file}\"\)\)\,\"dgamma\"\,list\(shape\=$mean_ind \* $mean_ind \/ \($sd_ind \* $sd_ind\)\, scale\=$sd_ind \* $sd_ind \/ $mean_ind\)\) | R --no-save | grep value | awk '{print $6}'`
echo "Ind: ", $mean_ind, $sd_ind , $pvalue_shapiro_ind, $pvalue_supremum_anderson_darling_ind
#echo sd_bootstrap, $sd_bootstrap
set mychi = `echo $file | awk -F "_" '{print $5 "/" substr($6,1,1)}'`
set sd_agg = `cat ${resname}${resnum}_choose5_w1530_p999_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//g ; s/ind_mutinfs.*mutinf\///g; s/mutinf.*ind_mutinf \=//; s/\n//g}' | grep "chi1\/chi2: $mychi" | awk '{print $12}' | tail -n 1`
echo "cat ${resname}${resnum}_choose5_w1530_p999_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//; s/ind_mutinfs.*m//}' | grep "chi1\/chi2: $mychi   utinf/ind_mutinf" | awk '{print $6 +$7}' | tail -n 1"
echo "cat ${resname}${resnum}_choose5_w1530_p999_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\\n//g ; s/ind_mutinfs.*mutinf\///g; s/mutinf.*ind_mutinf \=//; s/\\n//g}'"
set mean_agg = `cat ${resname}${resnum}_choose5_w1530_p999_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//g ; s/ind_mutinfs.*mutinf\///g; s/mutinf.*ind_mutinf \=//}' | grep "chi1\/chi2: $mychi" | awk '{print $9 + $10}' | tail -n 1`
set mean_agg_exc = `cat ${resname}${resnum}_choose5_w1530_p999_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//g; s/ind_mutinfs.*mutinf\///g; s/mutinf.*ind_mutinf \=//}' | grep "chi1\/chi2: $mychi " | awk '{print $9}' | tail -n 1`
echo mychi, $mychi
echo sd_agg, $sd_agg
echo $file, $sd_bootstrap, $sd_agg, $pvalue_shapiro,$pvalue_supremum_anderson_darling," Ind: ",$mean_ind,$sd_ind,$pvalue_shapiro_ind,$pvalue_supremum_anderson_darling_ind
echo $file, $mean_bootstrap, $mean_agg, $sd_bootstrap, $sd_agg,$pvalue_shapiro,$pvalue_supremum_anderson_darling," Ind: ",$mean_ind,$sd_ind,$pvalue_shapiro_ind,$pvalue_supremum_anderson_darling_ind," Excess: ",$mean_bootstrap_exc, $mean_agg_exc, $sd_bootstrap_exc, $pvalue_shapiro_exc,$pvalue_supremum_anderson_darling_exc >> compare_sds_choose3_w1530_p99.txt
set pdfoutfile = `echo $file | sed -e 's/.txt/_hist.pdf/'`
rm $pdfoutfile
echo pdf\(\"${pdfoutfile}\"\)\;hist\(t\(as.matrix\(read.table\(\"${file}\"\)\)\)\,main\=\"M.I. for Subsamples of ${file}\"\,xlab\=\"Mutual Information \(nats\)\" \)\;dev.off\(\) | R --no-save | tail -n 2 | head -n 1
set pdfoutfile_ind = `echo $ind_output_file | sed -e 's/.txt/_hist.pdf/'`
rm $pdfoutfile_ind
echo pdf\(\"${pdfoutfile_ind}\"\)\;hist\(t\(as.matrix\(read.table\(\"${ind_output_file}\"\)\)\)\,main\=\"Ind. M.I. for ${ind_output_file}\"\,xlab\=\"Mutual Information \(nats\)\" \)\;dev.off\(\) | R --no-save | tail -n 2 | head -n 1

end

end

cat > entsums_choose3_w1530_p99.txt <<EOF
$ent_sum_11
$ent_sum_12
$ent_sum_13
$ent_sum_14
$ent_sum_15
$ent_sum_16
$ent_sum_17
$ent_sum_18
$ent_sum_19
$ent_sum_20
EOF

set pdfoutfile_totent = entsums_choose3_w1530_p99_hist.pdf
rm $pdfoutfile_totent
echo pdf\(\"${pdfoutfile_totent}\"\)\;hist\(t\(as.matrix\(read.table\(\"entsums_choose3_w1530_p99.txt\"\)\)\)\,main\=\"Total Entropies from Bootstrap Samples \"\,xlab\=\"Mutual Information \(nats\)\" \)\;dev.off\(\) | R --no-save | tail -n 2 | head -n 1
set mean_totent = `echo mean\(t\(as.matrix\(read.table\(\"entsums_choose3_w1530_p99.txt\"\)\)\)\) | R --no-save | tail -n 2 | head -n 1 | awk '{print $2}'`
set sd_totent = `echo sd\(t\(as.matrix\(read.table\(\"entsums_choose3_w1530_p99.txt\"\)\)\)\) | R --no-save | tail -n 2 | head -n 1 | awk '{print $2}'`
echo "Mean Bootstrap Entropy: $mean_totent"
echo "Stdev Bootstrap Entropy: $sd_totent"
echo

#END
