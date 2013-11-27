#!/usr/bin/perl

$pdbfile=$ARGV[0];       #structure file
$resbfacfile=$ARGV[1];   #file with bfactors 
$resbfacfile2=$ARGV[2];  #file with bfactors 
$mut_res=$ARGV[3];       #residue having mutation (use GLY00 as default)
$outputprefix=$ARGV[4];  #outputprefix
$pdbfile =~ s/.pdb//;
$resbfacfile =~ s/.pdb//;
$resbfacfile2 =~ s/.pdb//;

open(STRUCTFILE,"${pdbfile}.pdb") or die "cannot open structure file ${pdbfile}.pdb!\n";
open(RESBFAC,"${resbfacfile}") or die "cannot open bfac file ${resbfacfile}!\n";
open(RESBFAC2,"${resbfacfile2}") or die "cannot open bfac file ${resbfacfile2}!\n";
open(OUTFILE,">${outputprefix}_resbfac.pdb");
open(PYMOLFILE,">${outputprefix}_resbfac.pml");

print OUTFILE "MODEL 1\n";
 

print PYMOLFILE "from pymol import cmd \n" ;
print PYMOLFILE "run ~cmcclend/mutinf/color_b.py \n";
print PYMOLFILE "load ${outputprefix}_resbfac.pdb, system \n" ;
print PYMOLFILE "cmd.remove('(all) and hydro') \n";
print PYMOLFILE "sele b0, b < 0.00001 and b > -0.00001 \n";
print PYMOLFILE "cmd.show('sticks','((byres (system))&n;ca,c,n,o,h)') \n"; #mainchain
#print PYMOLFILE "cmd.show('sticks','((byres (system))&(!(n;c,o,h|(n. n&!r. pro)))) and not b0')\n"; #sidechain
#print PYMOLFILE "cmd.show('sticks','system and (not b0)')\n";
print PYMOLFILE "cmd.hide('sticks'    ,'name O') \n ";
print PYMOLFILE "cmd.hide('lines'    ,'name O') \n ";
#print PYMOLFILE "cmd.spectrum('b',selection=('system'),quiet=0) \n";
print PYMOLFILE "color_b(selection='b<0',gradient='bw',mode='hist') \n"; 
print PYMOLFILE "color_b(selection='b>0',gradient='wr',mode='hist') \n"; 
print PYMOLFILE "cmd.color(5278, 'b0') \n";
print PYMOLFILE "cmd.disable('b0')  \n";
print PYMOLFILE "cmd.bg_color('white' ) \n"; 

$pymolstring = <<'END';

stored.list=[]
stored.list1=[]
stored.list2=[]
stored.list3=[]
cmd.iterate("(name ca)","stored.list1.append((resi,resn,color))")
cmd.iterate("(name n)","stored.list2.append((resi,resn,color))")
cmd.iterate("(name c)","stored.list3.append((resi,resn,color))")
#print stored.list

#python
i = 0

cmd.set_bond('stick_color', 5278 ,'system and (name C,N)')
python

for resnum_name in stored.list1: 
     cmd.select('this_res','resi '+str(resnum_name[0]))
     print "this res: "+str(resnum_name[0])
     cmd.select('this_res_N','this_res and name N')
     cmd.select('this_res_CA','this_res and name CA')
     cmd.select('this_res_C','this_res and name C')

     try:
       cmd.select('this_res_N','this_res and name N')
       stored.temp_list=[]
       cmd.iterate('this_res_N','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name N,CA)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name N,CA)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name N,CA)')
       #print "coloring N-CA: "+str(stored.temp_list[0])
     except IndexError:
     	    pass
     try:
       cmd.select('this_res_C','this_res and name C')
       stored.temp_list=[]
       cmd.iterate('this_res_C','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CA,C)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CA,C)')
       
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CA,C)')
       #print "coloring CA-C: "+str(stored.temp_list[0])
     except IndexError:
     	    pass

     #cmd.set_bond('stick_color',stored.list2[i][2],'this_res and (name N,CA)')
     #cmd.set_bond('stick_color',stored.list3[i][2],'this_res and (name CA,C)')

     try:
       cmd.select('this_res_CB','this_res and name CB')
       stored.temp_list=[]
       cmd.iterate('this_res_CB','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CA,CB)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CA,CB)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CA,CB)')
       print "coloring CA-CB: "+str(stored.temp_list[0])
     except IndexError:
     	    pass
     try:	  
       cmd.select('this_res_SG','this_res and name SG')
       stored.temp_list=[]
       cmd.iterate('this_res_SG','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CB,SG)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CB,SG)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CB,SG)')
       print "coloring CB-SG: "+str(stored.temp_list[0])
     except IndexError:
      	    pass
     try:
       cmd.select('this_res_OG','this_res and name OG')
       stored.temp_list=[]
       cmd.iterate('this_res_OG','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CB,OG)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CB,OG)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CB,OG)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_CG','this_res and name CG')
       stored.temp_list=[]
       cmd.iterate('this_res_CG','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CB,CG)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CB,CG)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CB,CG)')
       print "coloring CB-CG: "+str(stored.temp_list[0])
       
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_CG1','this_res and name CG1')
       stored.temp_list=[]
       cmd.iterate('this_res_CG1','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CB,CG1)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CB,CG1)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CB,CG1)')
     except IndexError:
      	    pass
     try:
       cmd.select('this_res_CG2','this_res and name CG2')
       stored.temp_list=[]
       cmd.iterate('this_res_CG2','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CB,CG2)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CB,CG2)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CB,CG2)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_OD1','this_res and name OD1')
       stored.temp_list=[]
       cmd.iterate('this_res_OD1','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CG,OD1)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CG,OD1)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CG,OD1)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_OD2','this_res and name OD2')
       stored.temp_list=[]
       cmd.iterate('this_res_OD2','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CG,OD2)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CG,OD2)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CG,OD2)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_CD','this_res and name CD')
       stored.temp_list=[]
       cmd.iterate('this_res_CD','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CG,CD)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CG,CD)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CG,CD)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_SD','this_res and name SD')
       stored.temp_list=[]
       cmd.iterate('this_res_SD','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CG,SD)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CG,SD)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CG,SD)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_CE','this_res and name CE')
       stored.temp_list=[]
       cmd.iterate('this_res_CE','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CD,CE)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CD,CE)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CD,CE)')
     except IndexError:
      	    pass
     try:
       cmd.select('this_res_NE','this_res and name NE')
       stored.temp_list=[]
       cmd.iterate('this_res_NE','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CD,NE)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CD,NE)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CD,NE)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_OE1','this_res and name OE1')
       stored.temp_list=[]
       cmd.iterate('this_res_OE1','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CD,OE1)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CD,OE1)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CD,OE1)')
     except IndexError:
      	    pass

     try:
       cmd.select('this_res_OE2','this_res and name OE2')
       stored.temp_list=[]
       cmd.iterate('this_res_OE2','stored.temp_list.append(color)')
       cmd.set_bond('stick_color',stored.temp_list[0],'this_res and (name CD,OE2)')
       cmd.set_bond('line_color',stored.temp_list[0],'this_res and (name CD,OE2)')
       if(stored.temp_list[0] != 5278): #if not zero
          cmd.show('sticks','this_res and (name CD,OE2)')
     except IndexError:
      	    pass

     i += 1


python end

END

print PYMOLFILE $pymolstring;
print PYMOLFILE "save ${outputprefix}_resbfac.pse \n";

#print PYMOLFILE "from pymol import cmd \n" ;
#print PYMOLFILE "load ${outputprefix}_resbfac.pdb, system \n" ;
#print PYMOLFILE "cmd.remove('(all) and hydro') \n";
#print PYMOLFILE "cmd.show('ribbon'    ,'system') \n";
#print PYMOLFILE "cmd.spectrum('b',selection=('system'),quiet=0) \n";
#print PYMOLFILE "sele b0, b < 0.00001 and b > -0.00001 \n";
#print PYMOLFILE "cmd.color(5278, 'b0') \n";
#print PYMOLFILE "cmd.disable('b0')  \n";
#print PYMOLFILE "cmd.bg_color('white' ) \n" ;




#print  "from pymol import cmd \n" ;
#print  "load ${outputprefix}_resbfac.pdb, system \n" ;
#print  "preset.b_factor_putty('system') \n" ;
#print  "sele b0, b < 0.00001 \n" ;
#print  "cmd.color(5278, 'b0') \n" ;
#print  "cmd.disable('b0') \n"     ;
#print  "cmd.bg_color('white') \n" ;

$bfacline = <RESBFAC>;
$bfacline2 = <RESBFAC2>;

if(! ($bfacline =~ /[0-9]+/)) #then this denotes a header line
{
  $bfacline = <RESBFAC>;
  $bfacline2 = <RESBFAC2>;
}
while($resline = <STRUCTFILE>)
{
    while(!(($resline =~ /ATOM/) || ($resline =~ /HETATM/))) 
    {
	#print $resline;
	if(!(eof STRUCTFILE))
	{
	    $resline = <STRUCTFILE>;
	}
	else
	{
	    exit;
	}
    }
    
    $tag    =substr($resline,0,6);
    $atnum  =substr($resline,6,5); 
    $atname =substr($resline,12,4);
    $resname=substr($resline,17,3);
    $chn    =substr($resline,21,1);
    $resnum =substr($resline,22,4);
    $x      =substr($resline,30,8);
    $y      =substr($resline,38,8);
    $z      =substr($resline,46,8);
    $resline =~ s/\n//;
    $occ = 1;
    $bfac = 0;
    $resname =~ s/\s+//;
    $resnum =~ s/\s+//;
    print "Looking for ${resname}${resnum}\n";
    if((($tag =~ /ATOM/) || ($tag =~ /HETATM/ )) && $resnum >= 1)
    {
	print $bfacline;
	while(!(eof RESBFAC) && !(($bfacline =~ /${resname}${resnum}${chn}/) || ($bfacline =~ /${resname}${resnum}/) || ($mut_res =~ /${resnum}/))) 
	{ 
	    last if(($resname =~ /HIE/) && ($bfacline =~ /HIS${resnum}/));
	    last if(($resname =~ /HID/) && ($bfacline =~ /HIS${resnum}/));
	    last if(($resname =~ /HIS/) && ($bfacline =~ /HIE${resnum}/));
	    last if(($resname =~ /HIP/) && ($bfacline =~ /HIS${resnum}/));
	    last if(($resname =~ /SEP/) && ($bfacline =~ /SER${resnum}/));
            last if(($resname =~ /TYS/) && ($bfacline =~ /TYR${resnum}/));
            #last if ($bfacline =~ /^${mut_res}/);
	    $bfacline = <RESBFAC>;
	    $bfacline2 = <RESBFAC2>;
	    print $bfacline;
	    print $bfacline2;
	} #loop
	$bfacline =~ s/KLDIV //; 
	$bfacline2 =~ s/KLDIV //; 
	($namenum, $bfac1, $bfac2, $bfac3, $bfac4, $bfac5, $bfac6) = split /\s+/, $bfacline;
	($namenum, $bfac1b, $bfac2b, $bfac3b, $bfac4b, $bfac5b, $bfac6b) = split /\s+/, $bfacline2;
	$bfac = 0;
	print "$bfac1 $bfac2 $bfac3 $bfac4 $bfac5 $bfac6\n";
	print "$bfac1b $bfac2b $bfac3b $bfac4b $bfac5b $bfac6b \n";
	
	$atname_underscores = $atname;
	$atname_underscores =~ s/ /_/g;
	print "atom name: ",$atname_underscores;
	#phi / psi
	if($atname =~ " N " ) {
	    $bfac = $bfac1b - $bfac1; #Phi
	}
	if($atname =~ " C " ) {
	    $bfac = $bfac2b - $bfac2; #Phi
	}

	#chi1
	if (! ($resname =~ "GLY" || $resname =~ "ALA")) {
	    if($atname =~ " CA " || $atname =~ " CB ") {
		$bfac = $bfac3b - $bfac3; #chi1
	    }
	}

	#chi2
	elsif ($resname =~ "SER" || $resname =~ "THR") {
	    if($atname =~ " OG ") {
		$bfac = $bfac4b - $bfac4; #chi2
	    }
	}
	elsif ($resname =~ "CYS" || $resname =~ "CYX") {
	    if($atname =~ " SG ") {
		$bfac = $bfac4b -$bfac4; #chi2
	    }
	}
	elsif (! ($resname =~ "VAL")) {
	    if($atname =~ " CG ") {
		$bfac = $bfac4b - $bfac4; #chi2
	    }
	}

	#chi3 
	if ($resname =~ "MET") {
	    if($atname =~ " SD ") {
		$bfac = $bfac5b - $bfac5; #chi3
	    }
	}
	elsif ($resname =~ "GLU" || $resname =~ "GLN" || $resname =~ "LYS" || $resname =~ "ARG") {
	    if($atname =~ " CD ") {
		$bfac = $bfac5b - $bfac5; #chi3
	    }
	}
	
	#chi4
	if ($resname =~ "LYS") {
	    if($atname =~ " CE ") {
		$bfac = $bfac6b - $bfac6; #chi4
	    }
	}
	elsif ($resname =~ "ARG") {
	    if($atname =~ " NE ") {
		$bfac = $bfac6b - $bfac6; #chi4
	    }
	}
	
	
	
	

	$bfac *= 100.0 ;
	
	if($bfacline =~ /^${mut_res}/)
	{
	    #$bfac = 0.0;
	}
	write OUTFILE;
    }
    
    
}

print OUTFILE "ENDMDL\n";

close(OUTFILE);
close(RESBFAC);
close(STRUCTFILE);




close(PYMOLFILE);
        

format OUTFILE =
@<<<<<@>>>>  @<< @>> @@>>>     @##.### @##.### @##.### @#.##@###.##
$tag,$atnum,$atname,$resname,$chn,$resnum,$x,$y,$z,$occ,$bfac
.

#END
