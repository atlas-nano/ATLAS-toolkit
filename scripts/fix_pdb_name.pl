#!/usr/bin/perl -w

use strict;

if (! @ARGV or $#ARGV < 2) {
    die "usage: fix_pdb_name.pl filebase startnum endnum\n";
}

my ($filebase, $startnum, $endnum) = @ARGV;

my ($curr_file);

if (! $startnum =~ /^\d+/) {
    die "Invalid starting number. Expected integer\n";
}

if (! $endnum =~ /^\d+/) {
    die "Invalid ending number. Expected integer\n";
}

my ($i, $j);

if ($startnum > $endnum) {
    my ($tempval) = $startnum;
    $startnum = $endnum;
    $endnum = $tempval;
}

for $i ($startnum .. $endnum) {
    my (@outarray) = ();
    $curr_file = $ARGV[0] . "." . $i;
    if (open INFILE, $curr_file) {
	print "Processing $curr_file....";
	while (<INFILE>) {
	    chomp;
	    push @outarray, $_;
	}
	close INFILE;

	my ($out_file) = $ARGV[0];
	if ($out_file =~ /(.+)_pdb/) {
	    $out_file = $1;
	}
	print "Done\n";
	$out_file .= "_" . $i . ".pdb";
	if (open OUTFILE, "> $out_file") {
	    print "Creating $out_file...";
	    for $j (0 .. $#outarray) {
		print OUTFILE "$outarray[$j]\n";
	    }
	    close OUTFILE;
	    print "Done\n";
	    system "/ul/maiti/src/util/scripts/fixcurvepdb.pl $out_file";  
	    system "rm -f $curr_file";
	}
    }
}
