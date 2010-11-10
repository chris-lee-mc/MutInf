#!/usr/bin/perl
#
#usage: WW_arrange_chis.pl file.dat
# where file.dat is a file with torsion angles as output by VMD such as 
# in ~mcclendon/ww/WWdShare/apo_1i6c/side_chains_dihedrals/*.dat
#
#

%long = (
          "A" => "ALA",
          "C" => "CYS",
          "D" => "ASP", 
          "E" => "GLU",
          "F" => "PHE",
          "G" => "GLY",
          "H" => "HIS",
          "I" => "ILE",
          "K" => "LYS",
          "L" => "LEU",
          "M" => "MET",
          "N" => "ASN",
          "P" => "PRO",
          "Q" => "GLN",
          "R" => "ARG",
          "S" => "SER",
          "T" => "THR",
          "V" => "VAL",
          "W" => "TRP",
          "Y" => "TYR" );

%NumChis = (
     "ALA"=>0, 
     "CYS"=>1, 
     "CYN"=>1,
     "CYX"=>1, 
     "ASP"=>2, 
     "AS4"=>2, 
     "GLU"=>3, 
     "GL4"=>3, 
     "GLH"=>3, 
     "PHE"=>2, 
     "GLY"=>0, 
     "HIS"=>2,
     "HIP"=>2, 
     "HIE"=>2, 
     "HID"=>2, 
     "ILE"=>2, 
     "LYS"=>4, 
     "LYP"=>4, 
     "LEU"=>2, 
     "MET"=>3, 
     "ASN"=>2, 
     "GLN"=>3,
     "PRO"=>1, 
     "ARG"=>4, 
     "SER"=>1, 
     "THR"=>1, 
     "VAL"=>1,
     "TRP"=>2, 
     "TYR"=>2, 
     "CTH"=>2  );


$file = $ARGV[0];

open(INFILE, "$file");

@stuff = split /_/, $file;

$resnum = $stuff[0];
$resname = substr($stuff[1],0,3);

print $stuff[0] . "      " . $stuff[1] . "\n";

#print $NumChis{"ALA"} . "\n";
print $resname . "     " . $NumChis{$resname} . "\n";

$num_chis_thisres = $NumChis{$resname};

#$torname[1] = "phi";
#$torname[2] = "psi";


for($run = 1; $run <=6; $run++) {
    for($i = 1; $i <= $num_chis_thisres; $i++ ) {
	$torname = "chi" . "$i";
	open($chis[$run][$i],">", "../run${run}/${torname}${resname}${resnum}.xvg");
    }
}
#if($num_chis_thisres >= 1) open($chis[1], "../run1/chi1${resname}${resnum}.xvg");
#if($num_chis_thisres >= 2) open($chis[2], "../run1/chi2${resname}${resnum}.xvg");
#if($num_chis_thisres >= 3) open($chis[3], "../run1/chi3${resname}${resnum}.xvg");
#if($num_chis_thisres >= 4) open($chis[4], "../run1/chi4${resname}${resnum}.xvg");


$, = "    ";


#SEPARTE INTO 6 BOOTSTRAP SAMPLES NEXT
$count = 0;
while(<INFILE>) {

    @myfields = split;
    $count++;
    $run = int($count / 154873.0) + 1;
    #print int($run) . "\n";
    for($i = 1; $i <= $num_chis_thisres; $i++ ) {

	print {$chis[$run][$i]}  $myfields[0], $myfields[$i] . "\n";

    }   
}

for($run = 1; $run <=6; $run++) {
    for($i = 1; $i <= $num_chis_thisres; $i++ ) {
	close $chis[$run][$i];
    }
}

#END
