#!/usr/bin/perl -w
use strict;
sub ProcessFile(@);

if (!@ARGV) {
    die "usage: $0 carnal_file [savefile] [starttime]\n";
}

my ($infile, $outfile, $start) = @ARGV;
my (@outarray, $i, $j, $numvals, $counter);

if (! $outfile) {
    $outfile = "carnal_rms.dat";
}

if (!$start) {
    $start = 0;
}

$numvals = 0;

ProcessFile($infile);

open OUTFILE, "> $outfile" or die "Cannot create $outfile: $!\n";

for $j (0 .. $numvals) {
    $counter = $start;
    for $i (0 .. $#outarray) {
	if ($outarray[$i][$j]) {
	    print OUTFILE "" . sprintf("%6d",$counter) .  "$outarray[$i][$j]\n";
	}
	$counter++;
    }
    print OUTFILE "\n";
}
close OUTFILE;

sub ProcessFile(@) {
    my ($fle) = $_[0];
    my (@datavals, $outstring, $i, $counter);

    -e $fle or die "Cannot locate $fle: $!\n";
    open INFILE, $fle or die "Cannot Open $fle: $!\n";
    $counter = 0;
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^L\d+\s+(.+)$/) {
	    @datavals = split /\s+/, $1;
	    if ($#datavals > 0) {
		if ($numvals < $#datavals) {
		    $numvals = $#datavals;
		}

		$outstring = "";
		for $i (0 .. $#datavals) {
		    if ($datavals[$i] =~ /\d+\.\d+/) {
			$outstring = sprintf("%10.7f", $datavals[$i]);
			$outarray[$counter][$i] = $outstring;
		    }
		}
		$counter++;
	    }
	}
    }
    close INFILE;

    if ($#outarray < 1) {
	die "ERROR: $fle is invalid\n";
    }
}
