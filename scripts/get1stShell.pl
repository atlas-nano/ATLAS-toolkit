#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);

sub init;
sub get1stShell;
sub printResults;

my ($FILES, $DATA);

$|++;
$FILES = &init;
&get1stShell($FILES);
&printResults($FILES);

sub printResults {
    my ($FILELIST) = $_[0];
    my ($i, $RESULTS, $num);

    for $i (@{ $FILELIST }) {
	if ($i->{FILE} =~ /(\d+\.\d+)/) {
	    $num = $1;
	}
	$RESULTS->{$num} = sprintf("%-8.3f %8.2f %12.6G # $i->{FILE}\n", $num, $i->{RADII}, $i->{FREQUENCY});
    }

    for $i (sort {$a<=>$b} keys %{ $RESULTS }) {
	print $RESULTS->{$i};
    }
}

sub get1stShell {
    my ($FILELIST) = $_[0];
    my ($i, $max, $numWSpace);

    for $i (@{ $FILELIST }) {
	print "Analyzing $i->{FILE}...\r";
	open INFILE, $i->{FILE} or die "ERROR: Cannot open $i->{FILE}: $!\n";
	$max = 0;
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /^(\d+\.\d+)\s+(\d+\.?\d*e?\-?\d*)/) {
		if ($2 > $max) {
		    $i->{RADII} = $1;
		    $i->{FREQUENCY} = $2;
		    $max = $2;
		}
	    }
	}
	close INFILE;
	die "ERROR: Data file $i is invalid\n" if (! defined($i->{RADII}));
	$numWSpace = length("Analyzing $i->{FILE}..");
    }
    printf "%-${numWSpace}s\n", "Analyzing files...Done";
}

sub init {
    my (%OPTS, @FILELIST, $rec, $findCmd);

    getopt('l',\%OPTS);

    die "usage: $0 -l data file location\n" if (! defined($OPTS{l}));

    print "Initializing...";
    $findCmd = "find $OPTS{l} -name '*.rdf' -print";
    open FINDCMD, "$findCmd |" or die "ERROR: Cannot execute $findCmd: $!\n";
    while (<FINDCMD>) {
	chomp;
	$rec = ();
	$rec->{FILE} = $_;
	push @FILELIST, $rec;
    }
    close FINDCMD;

    die "ERROR: No .dat files found while searching $OPTS{l}\n" if (! @FILELIST);

    print "\n";
    return \@FILELIST;
}
