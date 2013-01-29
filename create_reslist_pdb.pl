#!/usr//bin/perl


$pdbfile=$ARGV[0];       #structure file
$pdbfile =~ s/.pdb//;

open(STRUCTFILE,"${pdbfile}.pdb") or die "cannot open structure file ${pdbfile}.pdb!\n";

$gchi_counter = 0;


open(RESFILE,">${pdbfile}_reslist.reslist") or die "cannot create new reslist file ${pdbfile}_reslist.reslist";
@all_files = grep {/phi.*.xvg/ } readdir(RUN1);
print "files in ${subdir}/run1/ :" . (@all_files + 0) . "\n";


    
while($resline = <STRUCTFILE>)
{
    while(!(($resline =~ /ATOM/))) 
    {
	if( $gchi_counter > 1)
        {
	    if($resline =~ /MODEL/)  
	    {
		exit;
	    }
	    if($resline =~ /END/)
	    {
		exit;
	    }
	}
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
    $chain  =substr($resline,21,1);
    $resnum =substr($resline,22,4);
    $x      =substr($resline,30,8);
    $y      =substr($resline,38,8);
    $z      =substr($resline,46,8);
    
    $resline =~ s/\n//;
    $occ = 1;
    $bfac = 0;
    $resname =~ s/\s+//;
    $resnum =~ s/\s+//;
    if($resname =~ /ACE|NMA/) { next; }
    if($resname =~ /HIE/ || $resname =~ /HIP/ || $resname =~ /HID/) { $resname = "HIS"; }
    if($resname =~ /F3G/) { $resname = "CYS";}
    if($resname =~ /SEP/) { $resname = "SER";}
    if($resname =~ /TYS/) { $resname = "TYR";}
    if($resline =~ /^TER/) { next; }
    
	
    if(!($atname =~ / CA /)) 
    { 
	next;
    }


    $gchi_counter++;
    print("counter: $gchi_counter\n");
    print "${gchi_counter} ${resname} ${resnum}${chain}\n";
    print RESFILE "${gchi_counter} ${resname} ${resnum}${chain}\n";
       
    
    
    

}
close RESFILE;

#END
