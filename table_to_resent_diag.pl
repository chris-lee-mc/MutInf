#!/usr/bin/perl


%short = (
          "ALA" => "A",
          "CYS" => "C",
          "CYN" => "C",
          "CYX" => "C",
          "ASP" => "D",
          "GLU" => "E",
          "PHE" => "F",
          "GLY" => "G",
          "HIS" => "H",
          "HSD" => "H",
          "HSE" => "H",
          "HSP" => "H",
          "HIP" => "H",
          "ILE" => "I",
          "LYS" => "K",
          "LYP" => "K", 
          "LEU" => "L",
          "MET" => "M",
          "ASN" => "N",
          "PRO" => "P",
          "GLN" => "Q",
          "ARG" => "R",
          "SER" => "S",
          "THR" => "T",
          "VAL" => "V",
          "TRP" => "W",
	  "TYR" => "Y",
          "NAL" => "A",
          "CTH" => "T");


open(MYTABLE,'<',"${ARGV[0]}.txt");

open(MYDIAG,'>',"$ARGV[0]_resent.txt");

#read table

#find # residues
$numlines = 0;
$ent_tot=0;
while ($line = <MYTABLE>)
{
    $numlines++;
    if($numlines > 1)
    {
	$line =~ /([A-Z]+)([0-9]+)([A-Z]*) /;
	$res_name[$numlines - 1] = $1;
	$res_num[$numlines - 1] = $2;
	$chain[$numlines - 1] = $3;
	#print "$line" ;
	
	@values = split(/\s+/,$line);
	shift @values;
	#print @values;
	$index = 1;
	$res_ent[$numlines - 1] = 0;

	foreach $myvalue (@values)
	{

	    if($index == $numlines - 1)
	    {
		$res_ent[$numlines - 1] += $myvalue;
		print $myvalue;
	    }
	    else
	    {
		$res_ent[$numlines - 1] -= 0.0; # the other half belongs to the other residue
	    }
	    $index++;

	    #print "$myvalue ";
	    
	}
	$ent_tot += $res_ent[$numlines-1] ;
	#NEED TO FIX THE .RNA FILE ...SOMETHING IS NOT RIGHT, IT's SPLITTING UP DECIMALS AS WELL

    }
    print "\n";
}
$numres = $numlines - 1;
print "\n";
print "tot_ent: " + $ent_tot + "\n";

#write first row

print MYDIAG "RES\tRESENT\n"; # for .mrna format for cytoscape expression matrix

for($j=1;$j<=$numres;$j++)
{
    print MYDIAG "$res_name[$j]" . "$res_num[$j]" . "$chain[$j]" . "\t" . "$res_ent[$j]\n";
    
}
close MYDIAG;
