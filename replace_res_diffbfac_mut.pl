#!/usr/bin/perl

$pdbfile=$ARGV[0];       #structure file
$resbfacfile=$ARGV[1];   #bfactor file 
$resbfacfile2=$ARGV[2];  #bfactor file 
$mut_res=$ARGV[3];       #residue having mutation (use GLY00 as default)
$pdbfile =~ s/.pdb//;
$resbfacfile =~ s/.pdb//;

open(STRUCTFILE,"${pdbfile}.pdb") or die "cannot open structure file ${pdbfile}.pdb!\n";
open(RESBFAC,"${resbfacfile}") or die "cannot open bfac file ${resbfacfile}!\n";
open(RESBFAC2,"${resbfacfile2}") or die "cannot open bfac file ${resbfacfile}!\n";
open(OUTFILE,">${pdbfile}_resdiffbfac.pdb");

print OUTFILE "MODEL 1\n";
$bfacline = <RESBFAC>;
$bfacline2 = <RESBFAC2>;
if(! ($bfacline =~ /[0-9]+/)) #then this denotes a header line
{
  $bfacline = <RESBFAC>;
}
if(! ($bfacline2 =~ /[0-9]+/))
{
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
	while(!(eof RESBFAC) && !(($bfacline =~ /${resname}${resnum}/) || ($mut_res =~ /${resnum}/))) 
	{ 
	    last if(($resname =~ /HIE/) && ($bfacline =~ /HIS${resnum}/));
	    last if(($resname =~ /HID/) && ($bfacline =~ /HIS${resnum}/));
	    last if(($resname =~ /HIS/) && ($bfacline =~ /HIE${resnum}/));
	    #last if ($bfacline =~ /^${mut_res}/);
	    $bfacline = <RESBFAC>;
	    $bfacline2 = <RESBFAC2>;
	    print $bfacline;
	} #loop
	$bfacline =~ s/KLDIV //; 
	($namenum, $bfac) = split /\s+/, $bfacline;
	($namenum2, $bfac2) = split /\s+/, $bfacline2;
	$bfac = $bfac2 - $bfac;
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

format OUTFILE =
@<<<<<@>>>>  @<< @>> @@>>>     @##.### @##.### @##.### @#.##@###.##
$tag,$atnum,$atname,$resname,$chn,$resnum,$x,$y,$z,$occ,$bfac
.

#END
