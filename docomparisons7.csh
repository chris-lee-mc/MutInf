#!/bin/tcsh

cd dihedral_compare

# /home/cmcclend/dihedrals_3FJQ_atp_Mg2_HIP87_1.0/ /home/cmcclend/1CMK_apo_sapiens_HIP87_anton_nowat_1.0/ 
# on kirkwood

#foreach dirname1 (3FJQ_twoMg/dihedrals_3FJQ_atp_twoMg_HIP87_1.0  3FJQ_atp_Mg2_HIP87/dihedrals_3FJQ_atp_Mg2_HIP87_1.0/  1CMK_atp_Mg2_mouse_HIP87/dihedrals_1CMK_atp_Mg2_HIP87 1CMK_apo/dihedrals_1CMK_apo_HIP87_1.0 )
#foreach dirname1 #(3FJQ_atp_Mg2_HIP87/dihedrals_3FJQ_atp_Mg2_HIP87_1.0/ 
foreach dirname1 ( 3FJQ_twoMg/dihedrals_3FJQ_atp_twoMg_HIP87_1.0/ )
  #foreach dirname2 ( 3FJQ_atp_Mg2_HIP87/dihedrals_3FJQ_atp_Mg2_HIP87_1.0/  1CMK_atp_Mg2_mouse_HIP87/dihedrals_1CMK_atp_Mg2_HIP87 1CMK_apo/dihedrals_1CMK_apo_HIP87_1.0 ) 
  foreach dirname2 (  3FJQ_atp_twoMg_wPKI/dihedrals_3FJQ_atp_twoMg_wPKI_HIP87_1.0/ ) 
    
    if ($dirname1 != $dirname2) then
      set myprefix1=`echo ${dirname1} | sed -e 's/\/.*//' `
      set myprefix2=`echo ${dirname2} | sed -e 's/\/.*//' `
      mkdir ${myprefix1}_${myprefix2}
      cd ${myprefix1}_${myprefix2}

      foreach residue (`cat ../../3FJQ.reslist | awk '{print $2 $3}'` | sed -e 's/ //')
        python /home/cmcclend/mutinf/dihedral_analyze_compare.py ../../${dirname1} ../../${dirname2} $residue -d 14 # >& nohup_compare_${dirname1}_${dirname2}_out.txt &
      end
      cd ..
    endif
    
  end
end
cd ..

#END
