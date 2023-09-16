#!/usr/bin/perl -w
use strict;

if (! @ARGV or $#ARGV < 1) {
    die "usage: $0 datafile1 datafile2 ....\n";
}

sub GetValues(@);

my ($out_string, %NOEs, $is_valid);
my (@datafiles) = @ARGV;

$is_valid = 0;
for (@datafiles) {
    if (-e $_) {
	GetValues($_);
    }
}

print "Creating avgs.txt....";
if ($is_valid) {
    open OUTFILE, "> avgs.txt" or die "Cannot create avgs.txt: $!\n";
    for (keys %NOEs) {
	printf OUTFILE "%-30s", $_;
	print OUTFILE $NOEs{$_};
	print OUTFILE "\n";
    }
    close OUTFILE;
    print "Done\n";
}else {
    print "ERROR: No valid data found\n";
}

sub GetValues(@) {
    my ($infile) = $_[0];

    print "Obtaining data from $infile...";
    if (open INFILE, $infile) {
	while (<INFILE>) {
	    if ($_ =~ /^\s+\d+\s+(\d+\.\d+)\s+\d+\.\d+\s+(\w+\s*\S+\s*\-\s*\w+\s*\S+)/) {
		$NOEs{$2} .= sprintf("%8.3f", $1);
		$is_valid = 1;
	    }
	}
	close INFILE;
	print "Done\n";
    }else {
	print "Error getting data from $infile\n";
    }
}
