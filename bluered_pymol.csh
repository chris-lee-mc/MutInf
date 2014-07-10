#!/bin/tcsh
set prefix=`echo ${1} | sed -e 's/.pdb//'`
cat > ${prefix}.pml <<EOF
run /home/cmcclend/mutinf/color_b.py 
load ${prefix}.pdb, system 
preset.pretty('system')
color_b(selection='b<0',gradient='bw',mode='hist')
color_b(selection='b>0',gradient='wr',mode='hist')
sele b0, b < 0.00001 and b > -0.00001 
cmd.color(5278, 'b0')
cmd.disable('b0')
cmd.bg_color('white')
EOF

pymol ${prefix}.pml

#END
