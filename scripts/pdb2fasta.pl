#!/usr/local/bin/perl -w

($#ARGV == 1 && ($ARGV[0] eq "-x" || $ARGV[0] eq "-X")) 
    or die "Usage: bgf2fasta.pl -[x|X] pdb-file.\n";

open (PDB, "<$ARGV[1]") or die "Can't open file $ARGV[1].\n";

%short = (
	  "ALA" => "A",
	  "CYS" => "C",
	  "CYX" => "C",
	  "ASP" => "D",
	  "GLU" => "E",
	  "PHE" => "F",
	  "GLY" => "G",
	  "HIS" => "H",
	  "HSD" => "H",
	  "HSE" => "H",
	  "HSP" => "H",
   	  "ILE" => "I",
	  "LYS" => "K",
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
	  "TYR" => "Y" );
	
$cur_res = -1;  
foreach $line (<PDB>) {
    if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
	$res_num = substr($line, 23, 3);
	if ($res_num != $cur_res) {
	    $res_name = substr($line, 17, 3);
	    if (!$res_name) {          
		$res_name = "UNK";        #in case no res name
	    }
	    push @res_num_list, ($res_num);
	    push @res_name_list, ($res_name);
	    $cur_res = $res_num;
	}
	
    }
}

#$ARGV[1] =~ s/.pdb//;
#print ">$ARGV[1]", "\n";
for ($i = 0; $i <= $#res_name_list; $i++) {
    $a = $res_name_list[$i];
    if (!exists($short{$a})) {
	print "\n"; die "Unknown residue $a $res_num_list[$i]. No sequence output beyond this residue.\n";
    }
    if ($ARGV[0] eq "-X") {
	print $short{$a};
    }
    else {
	print lc $short{$a};
    }
    if (($i+1)%80 == 0) {
	print "\n";
    }
}

print "\n";	


