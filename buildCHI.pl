#!/usr/bin/perl
#1st argument is pdb file
#2nd argument is trajectory
open(A,"$ARGV[0]");
open(B,">","ptraj_tor.txt");
print B "trajin $ARGV[1]\n";
`mkdir run1`;
$firstres = -1;
  while (<A>) {
    if ($_ =~ /^ATOM/ && (($_ =~ /(C|O|S|N)G /) ||  ($_ =~ /(C|O|S|N)G1/  )))  {  
        $g = substr ($_,13,3);
        $g =~ s/ //g;
        $GG = 1;
    }
    elsif ($_ =~ /^ATOM/ && (($_ =~ /(N|S|C|O)D /) || ($_ =~ /(N|C|S|O)D1/)))  {  
        $d = substr ($_,13,3);
        $d =~ s/ //g;
        $DD = 1;
    }
    elsif ($_ =~ /^ATOM/ && (($_ =~ /(C|N|S|O)E /) || ($_ =~ /(C|N|S|O)E1/))) {
      $e = substr ($_,13,3);
      $e =~ s/ //g;
      $EE = 1;
    } 
    elsif ($_ =~ /^ATOM/ && (($_ =~ /(C|N|S|O)Z /) || ($_ =~ /(C|N|S|O)Z1/))) {
      $z = substr ($_,13,3);
      $z =~ s/ //g;
      $ZZ = 1;
    }
    elsif ($_ =~ /^ATOM/ && $_ =~ / O /)  {
    	$resnum = substr ($_,23,3);
        $resnum = $resnum*1;       
        $restype =  substr ($_,17,3);
	$i++;
        $n = ":".$resnum."\@N";
        $ca = ":".$resnum."\@CA";
        $cb = ":".$resnum."\@CB";
      
        if ($GG == 1) {
          $cg = ":"."$resnum"."\@".$g;
          $out1 = "run1/chi1".$restype.$resnum.".xvg";
          $name1 = $i."chi1";
          print B  "dihedral  $name1 $n $ca $cb $cg  out $out1 \n";         
          $GG = 0;
        }
         if ($DD == 1) {
          $cd = ":"."$resnum"."\@".$d;
	  $name2 = $i."chi2";
          $out2 = "run1\/chi2".$restype.$resnum.".xvg";
          print B  "dihedral  $name2 $ca $cb $cg $cd  out $out2 \n";
          $DD = 0;
        }
         if ($EE == 1) {
         $ce = ":"."$resnum"."\@".$e;
         $name3 = $i."chi3";
         $out3 = "run1\/chi3".$restype.$resnum.".xvg";
         print B  "dihedral  $name3  $cb $cg $cd $ce  out $out3 \n";
	 $EE = 0;
       }   
        if ($ZZ == 1) {
         $cz = ":"."$resnum"."\@".$z;
         $name4 = $i."chi4";
         $out4 = "run1\/chi4".$restype.$resnum.".xvg";
         print B  "dihedral  $name4  $cg $cd $ce $cz  out $out4 \n";
         $ZZ = 0;
       }
    
    
   }
   if ($_ =~ /^ATOM/ && $_ =~ / C /)  {
        $resnum = substr ($_,23,3);
	    
        $resnum = $resnum*1;  
	$plus =  $resnum +1;
	$minus = $resnum - 1;
        $restype =  substr ($_,17,3);

	$i++;
	$c = ":".$minus."\@C";
        $n = ":".$resnum."\@N";
        $ca = ":".$resnum."\@CA";
        $c2 = ":".$resnum."\@C";

        if ($firstres == -1) {
	    $firstres = 0;
	    $c = ":".$resnum."\@H1";
	}
        $name6 = $i."phi";
	$out6 = "run1\/phi".$restype.$resnum.".xvg";
	print B  "dihedral  $name6  $c $n $ca $c2  out $out6 \n";	

	$i++;
        $n = ":".$resnum."\@N";
        $ca = ":".$resnum."\@CA";
	$c = ":".$resnum."\@C";
        $n2 = ":".$plus."\@N";
        $name5 = $i."psi";
	$out5 = "run1\/psi".$restype.$resnum.".xvg";
	print B  "dihedral  $name5  $n $ca $c $n2  out $out5 \n";
	
	
   }
  } 

close(A);
close(B);
#now fix up last psi 
open(C,"ptraj_tor.txt");
$done = 0;
while(<C>)
{
  if(eof(C) == 1)
  {
      @myfields = split ;
      $myfields[4] =~ /[:]([0-9]+)[@]([A-Z])+/;
      $lastnum = $1;
  }
}
close(C);
print "$lastnum\n";
open(D,"ptraj_tor.txt");
print "a.prmtop\n";
$, = " ";
while(<D>)
{

    @myfields = split ;
    $myfields[5] =~ /[:]([0-9]+)[@]([A-Z,0-9]+)/;
    $mynum = $1;
    $mytext = "${lastnum}" . '@' . "OXT";
    if(int($mynum) == int($lastnum) + 1)
    {
	$myfields[5] =~ s/[:]([0-9]+)[@]([A-Z,0-9]+)/:${mytext}/ ;
	
    }

    print @myfields;
    print "\n";

}

