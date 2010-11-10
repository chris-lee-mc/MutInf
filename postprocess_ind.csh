#!/bin/tcsh

#cat ${1}_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//; }' | grep "mutinf/ind_mutinf" | awk '{print $2}' | sort -u 
set mychis = `cat ${1}_out.txt | sed -e '/chi1\/chi2/ { N; N; N; N; s/\n//g }' | grep "mutinf/ind_mutinf" | awk '{print $2}' | sort -u | sed -e 's/\//\\\//'`
echo $mychis
set n = 0
echo $n
foreach chi (`echo $mychis`) 
set n = `perl -e "print $n+1"`
echo $n
set chifield = `echo $chi | sed -e 's/\\\//_/'`
less ${1}_out.txt | egrep "ind_mutinfs|ind mutinfs" | sed -e 's/ind_mutinfs//' | sed -e 's/ind mutinfs//' | head -n $n | tail -n 1 > ${1}_${chifield}_ind.txt
end

#END

