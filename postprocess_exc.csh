#!/bin/tcsh

set mychis = `cat ${1}_out.txt | sed -e '/chi1\/chi2/ { N; s/\n//; p}' | grep "mutinf/ind_mutinf" | awk '{print $2}' | sort -u | sed -e 's/\//\\\//'`
echo $mychis
foreach chi (`echo $mychis`) 
set chifield = `echo $chi | sed -e 's/\\\//_/'`
less ${1}_out.txt | sed -e '/chi1\/chi2/ { N; s/\n//}' | grep "chi1\/chi2: $chi   mutinf/ind_mutinf" | awk '{print $5}' > ${1}_${chifield}_exc.txt
end

