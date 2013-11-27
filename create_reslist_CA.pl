#!/usr//bin/perl


$pdbfile=$ARGV[0];       #structure file
$subdir =$ARGV[1];       #directory underneath the current one in which run1 resides
$pdbfile =~ s/.pdb//;

open(STRUCTFILE,"${pdbfile}.pdb") or die "cannot open structure file ${pdbfile}.pdb!\n";

$gchi_counter = -2;

opendir(RUN1, "${subdir}/run1/");
open(RESFILE,">${pdbfile}_reslist.reslist") or die "cannot create new reslist file ${pdbfile}_reslist.reslist";
@all_files = grep {/phi.*.xvg/ } readdir(RUN1);
print "files in ${subdir}/run1/ :" . (@all_files + 0) . "\n";

closedir(RUN1);
opendir(RUN1, "${subdir}/run1/");
#print readdir(RUN1);
closedir(RUN1);
    
while($resline = <STRUCTFILE>)
{
    while(!(($resline =~ /ATOM/))) 
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
    $chain  =substr($resline,21,1);
    $resnum =substr($resline,22,4);
    $x      =substr($resline,30,8);
    $y      =substr($resline,38,8);
    $z      =substr($resline,46,8);
    if(!($atname =~ / CA /)) { next; }
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
	
    #opendir(RUN1, "${subdir}/run1/");
    @match_files = (); 
    @match_files2 = (); 
    @match_files3 = (); 
    @match_files4 = (); 
    #grep { /phi${resname}${gchi_counter}.xvg/ } readdir(RUN1);
    #closedir(RUN1);
    #if(!($match_files[0] =~ /phi${resname}${gchi_counter}.xvg/)) {  #try HIS fix first
    #   if($resname =~ /HIE/ || $resname =~ /HIP/ || $resname =~ /HID/) 
    #   { 
#	   $resname = "HIS"; 
    #	   if(!($match_files[0] =~ /phi${resname}${gchi_counter}.xvg/)) 
    #	   {  #try HIS fix first
#	   $resname = "HIE";
#	   }
#       }
#       
#    }
    if($resname =~ "HIS") 
    {
	$altname = "HIE";
	$altname2 = "HID";
	$altname3 = "HIP";
	print "looking also for HIE or HID\n";
    }
    else
    {
	$altname = $resname;
	$altname2 = $resname;
    }
#    $altname = $resname;
    while(!($match_files[0] =~ /phi${resname}${gchi_counter}.xvg/ || $match_files2[0] =~ /phi${altname}${gchi_counter}.xvg/ || $match_files3[0] =~ /phi${altname2}${gchi_counter}.xvg/ || $match_files4[0] =~ /phi${altname3}${gchi_counter}.xvg/ )) #while no files match
    {
	$gchi_counter++;
	print("counter: $gchi_counter\n");
	opendir(RUN1, "${subdir}/run1/");
	print "Looking for ${subdir}/run1/phi${resname}${gchi_counter}.xvg\n";
	@match_files = grep { /phi${resname}${gchi_counter}.xvg/  } readdir(RUN1);
	rewinddir(RUN1);
	@match_files2= grep { /phi${altname}${gchi_counter}.xvg/  } readdir(RUN1);
	rewinddir(RUN1);
	@match_files3= grep { /phi${altname2}${gchi_counter}.xvg/ } readdir(RUN1);
	rewinddir(RUN1);
	@match_files4= grep { /phi${altname3}${gchi_counter}.xvg/ } readdir(RUN1);
	#print $match_files3[0];
	closedir(RUN1);	
	print $match_files[0];
	print "\n";
	if($gchi_counter > (2*@all_files + 0)) 
	    {
		print "WARNING: couldn't find phi${resname}${gchi_counter}.xvg\n";
		last;
	    }
    }
    $chain =~ s/ //;
    print "${gchi_counter} ${resname} ${resnum}${chain}\n";
    print RESFILE "${gchi_counter} ${resname} ${resnum}\n";
    $g_chi_counter += 1;
    
}
close RESFILE;

#END
