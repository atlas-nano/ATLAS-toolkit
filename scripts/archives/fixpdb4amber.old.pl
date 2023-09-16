#! /usr/bin/perl -w
# fixpdb.pl - fixes the names of the bases in the namot2 pdb file
# to make it compatible with biograf and MSCFF4.2
# e.g changes GUA to G, THY to  T etc.

# Variable Declaration Section

$outputarry[0] = "";
$flenm = "";
$in_text  = "";
# Opens the file and gets the input

if (!@ARGV) {
    die "usage: fixpdb.pl pdbfile.pdb\n";
}

$flenm = $ARGV[0];

if (!open(PDBFILE, $flenm)) {
    die "Cannot open $flenm, $!\n";
} else {
    while (<PDBFILE>) {
	$in_text = $_;
	if ($in_text =~ /(GUA|THY|ADE|CYT)/) {
	    push @outputarry, $in_text;
	} else {
	    if ($in_text =~ /HE/) {
		push @outputarry, "TER\n";
	    }
	}
    }
    close PDBFILE;
}

if (!open (OUTFILE, "> " . $flenm)) {
    die "Cannot open $flenm";
} else {
    for ($i=1; $i<=$#outputarry; $i++) {
	print OUTFILE $outputarry[$i];
    }
}

close OUTFILE;

print "Done\n";
