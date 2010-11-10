#!/bin/tcsh

less ${1}_out.txt | grep "Entropy" | awk '{print $2}'  > ${1}_entropy.txt
less ${1}_out.txt | grep "Entropy" | awk '{print $4}'  > ${1}_sd_ent.txt

